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
#include "Constraints/GasConstraint.hpp"

using namespace std;

class GasInjector
{
public:
    GasInjector(glm::dvec2 pos, GasConstraint *constr, double particles_per_sec, bool isOpen) {
        position = pos;
        speed = particles_per_sec;
        time = 0;
        constraint = constr;
        boundary = isOpen;
    }

    std::vector<Particle*> *getParticles() {
        return &particles;
    }

    glm::dvec2 getPosition() { 
        return position; 
    }

    GasConstraint *getConstraint() {
        return constraint;
    }

    void advanceTime(std::vector<Particle*> *in_particles, double elapsed) {
        time += elapsed;

        // Adding a new gas/smoke particle every 1/speed second 

        // Step 1: Creating new particle(s) based on type of the boundary
        while(time >= 1.0/speed) {
            time -= 1.0/speed;

            // The mass of the newly injected particles will be much less compared to other fluids
            Particle *new_particle = new Particle(position[0], position[1], 0.1, glm::vec3(0.99f, 0.99f, 0.99f), GAS);

            particles.push_back(new_particle);

            // To generate plumes of smoke, we inject fluid particles and smoke particles together (Section 7.2.2)
            if(boundary) {
                new_particle = new Particle(position[0], position[1], 1, glm::vec3(0.8f, 0.5f, 0.0f), GAS);
                constraint->addParticle(new_particle, in_particles->size());
                in_particles->push_back(new_particle);
            }
        }

        // Step 2: Assigning velocity of the particles baded on Equation (28) of Muller's "Unified Particles Physics" Paper
        for(Particle *p: particles) {
            if(p->ph == GAS || p->ph == FLUID) {
                p->v = glm::dvec2();
                double acc = 0;

                for(Particle *p2: *in_particles) {

                    glm::dvec2 r = p->p - p2->p;
                    double poly6_term = poly6(glm::dot(r,r));

                    p->v += p2->v * poly6_term;
                    acc += poly6_term;
                }

                if(acc > 0)
                    p->p += p->v * elapsed / acc;
            }
        }
    }

    

    double time;                        // Used for the elapsed time interval (timestep)
    double speed;                       // Number of particles that will be emitted per second
    bool boundary;                      // true for open boundaries and false for closed ones
    glm::dvec2 position;                // Position of the injector
    std::vector<Particle*> particles;  
    GasConstraint *constraint;
    
};
