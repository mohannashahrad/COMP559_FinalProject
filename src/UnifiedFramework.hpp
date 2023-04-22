/**
 * UnifiedFramework class 
 * @author Mohanna Shahrad 
 */

#include <stdlib.h>
#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec2.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include "GLSL.h"

#include "Particle.hpp"
#include "GasInjector.hpp"
#include "Constraints/FluidConstraint.hpp"

using namespace std;

// Number of solver iterations per time step 
// Tuned this and seems like 3 is enough. 
// More than 4 makes the visuals slower
#define SOLVER_ITRS 3

// Gravity scaling factor for gases
#define GRAVITY_REDUCTION_ALPHA -0.3
#define EPSILON 0.0001

// Scenes that are supported by the framework
enum Scene {
    GAS_OPEN,
    GAS_CLOSED,
    TWO_FLUID
};

enum Boundary {
    OPEN,
    CLOSED,
};

class UnifiedFramework
{
public:
    UnifiedFramework();
    virtual ~UnifiedFramework();

    // Initializers for test scenes
    void init(Scene type);
    void gasOpen();
    void gasClosed();
    void twoFluid();

    void advanceTime(double elapsed);
    void display();
    void resize(int width, int height);
    void showGrid();
    void unshowGrid();

    // Debug information and flags
    int getNumParticles();
    glm::dvec2 gravity;
    bool isGridShown;

private:

    // Reset the simulation
    void clear();

    // Functions for different types of particles
    GasConstraint *createGasConstraint(std::vector<Particle*> *particles, double density, Boundary boundary);
    FluidConstraint *createFluidConstraint(std::vector<Particle *> *in_particles, double density);
    void addGasInjector(glm::dvec2 posn, GasConstraint *gs, double particlesPerSec, bool isOpen);

    void displayGrid();
    void displayBoundaries();
    void displayParticles();
    void drawBodies();
    void displaySmoke();
    void fixContactConstraints();

    // Counts for iterative particle solver [parameter n in Alg 1 of the paper]
    int *sim_n;

    // Storage of global particles, rigid bodies, and general constraints
    std::vector<Particle*> particles;
    std::vector<GasInjector*> gasInjectors;
    std::map<ConstraintGroup, std::vector<Constraint *>> globalConstraints;

    // Drawing and boundary information
    glm::ivec2 dimensions;
    glm::dvec2 xBnd, yBnd;
};
