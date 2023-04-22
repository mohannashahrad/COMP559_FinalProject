/**<!-------------------------------------------------------------------->
   Implementation of poly6 and spikyGradient smoothing kernels expressed in 
   Muller's "Particle-Based Fluid Simulation for Interactive Applications" paper.

   Used for Fluids and Gas Constraints. 

   Poly6 is mainly used for density estimation while spikyGradient is used for gradient calculations.
   <!-------------------------------------------------------------------->**/
#ifndef SMOOTHINGKERNELS_H
#define SMOOTHINGKERNELS_H

// TODO: Play with values of H
#define H 2.
#define H2 4. 
#define H6 64.
#define H9 512.

// Look at equation (13) of Position-based fluids 
// The values are set to those suggested in the paper
#define N_TERM 4
#define K_TERM 0.1
#define DELTA_Q_TERM_FLUID 0.2
#define DELTA_Q_TERM_GAS 0.25

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/fwd.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/vec2.hpp>
#include <cmath>


double poly6(double r_squared);

// Equation (21) of Muller's "Particle-Based Fluid Simulation for Interactive Applications" paper
glm::dvec2 spikyGradient(const glm::dvec2 &r_vect, double r);

#endif // SMOOTHINGKERNELS_H
    