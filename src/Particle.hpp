#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <map>

#include <vector>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/vec2.hpp>
#include "SmoothingKernels.h"

typedef glm::vec<2, double> vec2;

using namespace std;

#define EPSILON 0.0001

// Size of the particles
#define PARTICLE_RAD 0.12
#define PARTICLE_DIAMETER 0.24

// Phase for particles
enum Phase {
    FLUID,
    GAS,
    RIGID,
    NUM_PHASES
};

// Groups of constraints
enum ConstraintGroup {
    CONTACT,
    GENERAL,
    NUM_CONSTRAINT_GROUPS
};

/**
 * Particle class that contains particle properties (e.g., mass), 
 * initial positions and velocities, current position and velocities 
 * and a force accumulator for computing the total force acting on this particle.
 * @author kry
 */
class Particle {
public:
    /** Identifies this particles position in the particle list */
    int index;

    bool pinned = false;

    glm::vec3 color{ 0.0f, 0.95f, 0.0f };

    double mass = 1;

    double i_mass;  // Inverse Mass (M^-1)
    double t_mass;  // Temporary Mass

    // Used in boundary constraitns
    double sFriction, kFriction;     // Friction coefficients (Static and Kinetic frictions)
    Phase ph;       // Phase of the particle (Look at the paper for more definition)  

    vec2 p;         // Position
    vec2 ep;        // Updated Position with velocity (an estimate value)
    vec2 v;
    vec2 vTemp;
    vec2 p0;
    vec2 v0;
    vec2 f;

    /** Default constructor */
    Particle() {}

    /**
     * Creates a particle with the given position, and velocities.
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
     Particle( double x, double y, double vx, double vy): vTemp(0,0), v(0,0) {
        p0 = vec2(x, y);
        v0 = vec2(vx, vy);
        reset();
    }

    /* Creates a particle with the given position, mass, and phase */
    Particle( double x, double y, double mass, Phase phase = RIGID): vTemp(0,0), v(0,0), ph(phase) {
        p0 = vec2(x, y);
        v0 = vec2(0,0);
        sFriction = 0;
        kFriction = 0;

        ep = glm::dvec2();
    
        if (mass <= 0) {
            i_mass = -mass;
        } else {
            i_mass = 1. / mass;
        }

        t_mass = i_mass;
        reset();
    }

    /* Creates a particle with the given position, mass, color, and phase */
    Particle( double x, double y, double mass, glm::vec3 newCol ,Phase phase = RIGID): vTemp(0,0), v(0,0), ph(phase), color(newCol) {
        
        p0 = vec2(x, y);
        v0 = vec2(0,0);
        sFriction = 0;
        kFriction = 0;

        ep = glm::dvec2();
    
        if (mass <= 0) {
            i_mass = -mass;
        } else {
            i_mass = 1. / mass;
        }

        t_mass = i_mass;
        reset();
    }

    /* Creates a particle with the given position, velocity, mass, and phase  */
    Particle( double x, double y, double vx, double vy, double mass, Phase phase = RIGID): vTemp(0,0), v(0,0), ph(phase) {
        p0 = vec2(x, y);
        v0 = vec2(vx, vy);

        // init values
        sFriction = 0;
        kFriction = 0;

        ep = glm::dvec2();
    
        if (mass <= 0) {
            i_mass = -mass;
        } else {
            i_mass = 1. / mass;
        }

        t_mass = i_mass;
        reset();
    }

    /**
     * Resets the position of this particle
     */
    void reset() {
        p = p0;
        v = v0;
        f = vec2(0, 0);
    }

    /**
     * Clears all forces acting on this particle
     */
    void clearForce() {
        f = vec2(0, 0);
    }

    /**
     * Adds the given force to this particle
     * @param force
     */
    void addForce(vec2 force) {
        f += force;
    }

    void setColor(glm::vec3 newCol) {
        color = newCol;
    }

    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    float distance(double x, double y) {
        vec2 diff = p - vec2(x, y);
        return (float) sqrt( diff.x*diff.x + diff.y*diff.y);
    }

    void update_pos_epsilon_check() {

        // Epsilon check
        if (glm::length(ep - p) < EPSILON) {
            v = vec2(0,0); 
            return;
        }

        // Update position
        p = ep;
    }

    // This is explained in Equation (21) of the paper
    // This step is only useful during the contact processing phase
    void scaleMass() {
        if (i_mass != 0.0) {
            t_mass = 1. / ((1. / i_mass) * exp(-p.y));
        } else {
            t_mass = 0.0;
        }
    }
};

// Super class of all constraint types
class Constraint {
public:
    Constraint() {}

    // Keeping the functions "virtual" as we want to re-define them in the subclasses
    virtual ~Constraint() {}

    // For iterative constraint solver
    virtual void project(std::vector<Particle*> *e_particles, int *counts) = 0;

    virtual void updateN(int *sim_n) = 0;
};
