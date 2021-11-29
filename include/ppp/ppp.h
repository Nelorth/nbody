#include <stdbool.h>
#include <stdlib.h>
#include <sys/time.h>

static const long double G = 6.674e-11;  // the gravitational constant

/**
 * Representation of bodies in memory.
 */
typedef struct {
    long double mass;    /* in kilograms */
    long double x, y;    /* position in meters */
    long double vx, vy;  /* velocity in meters per second */
} body;

/**
 * Container for simulation parameters.
 */
struct TaskInput {
    // the number of bodies
    int nBodies;

    // the bodies (also used for output)
    body* bodies;

    // whether use Newton's third law for local computations
    bool newton3local;

    // whether use Newton's third law globally
    bool newton3;

    // number of simulation steps to perform
    int nSteps;

    // the length of a time step of the simulation (in seconds)
    long double deltaT;

    // number of steps between two images
    // 0 means no image output
    int imageStep;

    // print some debug outputs during computation
    bool debug;
};

/**
 * Performs the nbody simulation sequentially in a single process.
 */
void compute_single(struct TaskInput* TI);

/**
 * Performs the nbody simulation in a distributed fashion using MPI and OpenMP.
 */
void compute_parallel(struct TaskInput* TI);

/**
 * Saves an image with ID <imgNum> depicting the given bodies in a 2D grid.
 */
void saveImage(int imgNum, const body* bodies, int nBodies);

/**
 * Returns the number of seconds since 1970-01-01T00:00:00.
 */
inline static double seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + ((double) tv.tv_usec) / 1000000;
}
