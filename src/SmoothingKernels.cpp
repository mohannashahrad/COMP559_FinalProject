#include "SmoothingKernels.h"

// Used for density estimation
double poly6(double r_squared) {
    if(r_squared >= H2) return 0;
    return (315.0 / (64.0 * M_PI * H9)) * ((H2 - r_squared) * (H2 - r_squared) * (H2 - r_squared));
}

// Equation (21) of Muller's "Particle-Based Fluid Simulation for Interactive Applications" paper
// Used for gradient Calculation
glm::dvec2 spikyGradient(const glm::dvec2 &r_vect, double r) {
    if(r >= H || r == 0) return glm::dvec2();

    // spiky = (15.0 / (M_PI * H6)) * ((H - r) * (H - r) * (H - r))
    return -glm::normalize(r_vect) * (45. / (M_PI * H6)) * (H - r) * (H - r);
}
    