#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#include "ppp/ppp.h"
#include "ppp_pnm/ppp_pnm.h"

/**
 * Computes the force body j exercises on body i.
 * The acceleration is returned in *ax and *ay.
 */
static void compute(body* bodies, int i, int j,
                    long double* ax, long double* ay) {
    long double aij_x, aij_y, dx, dy, r3;
    dx = bodies[j].x - bodies[i].x;
    dy = bodies[j].y - bodies[i].y;
    r3 = powl(sqrtl(dx * dx + dy * dy), 3);
    aij_x = dx / r3;
    aij_y = dy / r3;
    *ax = aij_x * bodies[j].mass;
    *ay = aij_y * bodies[j].mass;
}

/**
 * Updates (x,y) and (vx,vy) of body b from its
 * current position, velocity and acceleration.
 */
static void update(body* b, double deltaT, long double ax, long double ay) {
    long double dvx, dvy;
    dvx = ax * (G * deltaT);
    dvy = ay * (G * deltaT);
    b->x += (b->vx + dvx / 2) * deltaT;
    b->y += (b->vy + dvy / 2) * deltaT;
    b->vx += dvx;
    b->vy += dvy;
}

/**
 * Determines length and offset in a problem buffer of size <size>
 * for process <rank> where the total number of processes is <np>.
 */
static void portionize(int size, int np, int rank, int* length, int* offset) {
    int portion = size / np;
    int residue = size % np;
    if (rank < residue) {
        *length = portion + 1;
        *offset = (*length) * rank;
    } else {
        *length = portion;
        *offset = (*length) * rank + residue;
    }
    if (*length == 0) {
        *offset = 0;
    }
}

/**
 * Determines x and y accelerations exercised by body j on body i.
 */
void acc_foo(body* bodies, int i, int j, long double* accX, long double* accY) {
    long double ax, ay;
    compute(bodies, i, j, &ax, &ay);
    accX[i] += ax;
    accY[i] += ay;
}

/**
 * Determines x and y accelerations exercised by body j on body i and
 * uses those to symmetrically set accelerations of body i on body j.
 */
void acc_bar(body* bodies, int i, int j, long double* accX, long double* accY) {
    long double ax, ay;
    compute(bodies, i, j, &ax, &ay);
    accX[i] += ax;
    accY[i] += ay;
    long double ratio = bodies[i].mass / bodies[j].mass;
    accX[j] -= ratio * ax;
    accY[j] -= ratio * ay;
}

/**
 * Performs the nbody simulation in a distributed fashion using MPI and OpenMP.
 */
void compute_parallel(struct TaskInput* TI) {
    int np, self;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &self);

    const bool debug = TI->debug;
    const long double deltaT = TI->deltaT;
    const int nSteps = TI->nSteps;
    const int imageStep = TI->imageStep;
    const int nBodies = TI->nBodies;
    body* bodies = TI->bodies;

    MPI_Datatype body_type, pos_type;
    MPI_Type_contiguous(5, MPI_LONG_DOUBLE, &body_type);
    MPI_Type_commit(&body_type);
    MPI_Type_create_indexed_block(1, 2, (int[]){1}, MPI_LONG_DOUBLE, &pos_type);
    MPI_Type_create_resized(pos_type, 0, sizeof(body), &pos_type);
    MPI_Type_commit(&pos_type);

    int counts[np], displs[np];
    #pragma omp parallel for
    for (int p = 0; p < np; p++) {
        portionize(nBodies, np, p, &counts[p], &displs[p]);
    }
    int myoffset = displs[self];
    int mylength = counts[self];

    long double accX[nBodies], accY[nBodies];

    for (int step = 0; step < nSteps; step++) {
        // save an image snapshot every <imageStep> steps
        if (self == 0 && imageStep > 0 && step % imageStep == 0) {
            saveImage(step / imageStep, bodies, nBodies);
        }

        if (self == 0 && debug) {
            printf("%d\r", step);
        }

        // initialize this step's accelerations
        #pragma omp parallel for
        for (int k = 0; k < nBodies; k++) {
            accX[k] = accY[k] = 0;
        }

        if (TI->newton3) {
            // implementation with Newton's third law used globally

            /*
             * reduction on whole arrays turned out to be faster than
             * initiating parallelism every iteration and reducing individual
             * fields manually into temporary local variables.
             */
            #pragma omp parallel for reduction(+:accX,accY)
            for (int i = myoffset; i < myoffset + mylength / 2; i++) {
                for (int j = 0; j < i; j++) {
                    acc_bar(bodies, i, j, accX, accY);
                }
            }

            int p = (np - 1) - self;  // mirror process
            #pragma omp parallel for reduction(+:accX,accY)
            for (int i = displs[p]+counts[p]/2; i < displs[p]+counts[p]; i++) {
                for (int j = 0; j < i; j++) {
                    acc_bar(bodies, i, j, accX, accY);
                }
            }

            // sum accelerations together over all processes
            MPI_Request requests[2];
            MPI_Iallreduce(MPI_IN_PLACE, &accX, nBodies, MPI_LONG_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD, &requests[0]);
            MPI_Iallreduce(MPI_IN_PLACE, &accY, nBodies, MPI_LONG_DOUBLE,
                    MPI_SUM, MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
        } else if (TI->newton3local) {
            // implementation with Newton's third law for local computations
            #pragma omp parallel for reduction(+:accX,accY)
            for (int i = myoffset; i < myoffset + mylength; i++) {
                for (int j = 0; j < nBodies; j++) {
                    if (j >= myoffset && j < myoffset + mylength) {
                        if (j < i) {
                            acc_bar(bodies, i, j, accX, accY);
                        }
                    } else {
                        acc_foo(bodies, i, j, accX, accY);
                    }
                }
            }
        } else {
            // straightforward implementation without calculation savings
            #pragma omp parallel for
            for (int i = myoffset; i < myoffset + mylength; i++) {
                for (int j = 0; j < nBodies; j++) {
                    if (i == j) {
                        continue;
                    }
                    acc_foo(bodies, i, j, accX, accY);
                }
            }
        }

        // update this process' associated bodies
        #pragma omp parallel for
        for (int i = myoffset; i < myoffset + mylength; i++) {
            update(&bodies[i], deltaT, accX[i], accY[i]);
        }

        // sync body positions over all processes
        MPI_Allgatherv(MPI_IN_PLACE, 0, pos_type, bodies, counts, displs,
                       pos_type, MPI_COMM_WORLD);
    }

    // save a final snapshot if <imageStep> divides <nSteps>
    if (self == 0 && imageStep > 0 && nSteps % imageStep == 0) {
        saveImage(nSteps / imageStep, bodies, nBodies);
    }

    if (debug) {
        printf("\n");
    }

    // collect final result on root
    if (self == 0) {
        MPI_Gatherv(MPI_IN_PLACE, 0, body_type, bodies, counts, displs,
                    body_type, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv(bodies + myoffset, mylength, body_type, bodies, counts,
                    displs, body_type, 0, MPI_COMM_WORLD);
    }

    // clean up custom types
    MPI_Type_free(&body_type);
    MPI_Type_free(&pos_type);
}
