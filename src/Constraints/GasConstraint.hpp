// Fluid-solid coupling constant
// Check equation (27) of Muller's "Unified Particle Physics for Real-Time Applications" paper
#define S_SOLID_COUPLING 0.4

#include "../Particle.hpp"


// Basics of this class is very similar to FluidConstraint
class GasConstraint : public Constraint
{
public:
    GasConstraint(double density, std::vector<int> *particles, bool isOpen) : Constraint() {
        open = isOpen; 
        rest_density = density;
        num_Particles = particles->size();

        // Indices of particles
        neighbors = new std::vector<int>[num_Particles];

        delta_values = new glm::dvec2[num_Particles];

        for (int i = 0; i < num_Particles; i++) {
            particle_indices.push_back((*particles)[i]);
        }
    }
    
    virtual ~GasConstraint() {
        delete[] neighbors;
    }

    // Based on Algorithm 1 of Muller's Position Based Fluids
    void project(std::vector<Particle*> *particles, int *sim_n) {
        lambda_values.clear();
        
        // Compute Lambda Values (Equation 11 of Mullers's Position-based Fluids)
        for (int i = 0; i < particle_indices.size(); i++) {
            
            neighbors[i].clear();

            // Index of that particle
            int j = particle_indices[i];
            Particle *p_j = particles->at(j);

            // density estimate of p_j
            double dens_j = 0.0;
            double denominator = 0.0;

            // Finding neighbors
            for (int k = 0; k < particles->size(); k++) {

                if (k != j) {
                    Particle *p_k = particles->at(k);

                    // Skip fixed particles
                    if (p_k->i_mass == 0) continue;

                    glm::dvec2 r_vect = p_j->ep - p_k->ep;
                    double r = glm::dot(r_vect, r_vect);    // || r_vect || = r

                    if (r < H2) {  // This is a neighbor
                        neighbors[i].push_back(k);
                        double poly6_term = poly6(r);

                        if (p_k->ph == RIGID) {
                            dens_j += ((poly6_term / p_k->i_mass) * S_SOLID_COUPLING);
                        } else {
                            dens_j += (poly6_term / p_k->i_mass);
                        }

                        glm::dvec2 gradient_vect = calculateGrad(particles, i, k);
                        denominator += glm::dot(gradient_vect, gradient_vect);
                    }

                // If next particle is the particle with index i
                } else {
                    neighbors[i].push_back(k);
                    dens_j += poly6(0) / p_j->i_mass;
                }
            }

            glm::dvec2 gradient_vect = calculateGrad(particles, i, j);
            denominator += glm::dot(gradient_vect, gradient_vect);
            double dens_ratio = (dens_j/rest_density);

            // Enforcing of f_drag [Equation 29 of Paper]
            // Using k of -2
            if(open) {
                p_j->addForce(p_j->v * (1.0-dens_ratio) * -2.0);
            }

            double relaxation = 0.05;   // For gamma correction
            lambda_values[j] = (-(dens_ratio - 1.) / (denominator + relaxation));
        }

        // Compute Delta Values (Equation 12 of Mullers's Position-based Fluids)
        // NOTE: spikyGradient kernel is used here instead of poly6 as Position-based Fluids
        for (int i = 0; i < particle_indices.size(); i++) {
            // Value of delta would be summed in the inner loop and updated in the outer loop 
            glm::dvec2 delta = glm::dvec2();
            glm::dvec2 f_vortocity = glm::dvec2();

            int j = particle_indices[i];
            Particle *p_j = particles->at(j);

            // Iterating through all neighbors of p_j
            for (int k = 0; k < neighbors[i].size(); k++) {

                int n = neighbors[i][k];

                if (j != n) {
                    // p_n is the neighbor of p_j
                    Particle *p_n = particles->at(n);

                    // Implementing Equation (13-14) of Muller's PBF
                    glm::dvec2 r_vect = p_j->ep - p_n->ep;
                    double r = glm::length(r_vect);

                    glm::dvec2 spiky_grad_term = spikyGradient(r_vect, r);

                    double sCorr = -K_TERM * pow((poly6(r * r) / poly6(DELTA_Q_TERM_GAS * DELTA_Q_TERM_GAS * H * H)), N_TERM);
                    delta += (lambda_values[j] + lambda_values[n] + sCorr) * spiky_grad_term;

                    // Applying Vortocity
                    // Look at section 7.2.3 of the unified physics paper
                    glm::dvec2 w = spikyGradient(r_vect, glm::dot(r_vect,r_vect)) * p_n->v;
                    glm::dvec3 crossProduct = glm::cross(glm::dvec3(0,0,glm::length(w)), glm::dvec3(r_vect.x, r_vect.y, 0));
                    
                    f_vortocity += glm::dvec2(crossProduct.x, crossProduct.y) * poly6(glm::dot(r_vect,r_vect));
                }
            }

            delta_values[i] = (delta / rest_density);
            p_j->f += f_vortocity;
        }

        for (int i = 0; i < particle_indices.size(); i++) {
            int j = particle_indices[i];
            Particle *p_j = particles->at(j);
            p_j->ep += delta_values[i] / ((double) neighbors[i].size() + sim_n[j]);
        }
    }
    
    void addParticle(Particle *p, int index) {
        particle_indices.push_back(index);

        delete[] neighbors;
        delete[] delta_values;

        num_Particles++;

        // Reallocare neighbors and  delta_values
        neighbors = new std::vector<int>[num_Particles];
        delta_values = new glm::dvec2[num_Particles];
        
    }

    // For Gradient Calculation spikyGradient is used as explained in Mulller's PBF
    glm::dvec2 calculateGrad(std::vector<Particle *> *particles, int idx1, int idx2){

        Particle *p1 = particles->at(particle_indices[idx1]);
        Particle *p2 = particles->at(idx2);

        // difference between p1 and p2
        glm::dvec2 r_vect = p1->ep - p2->ep;
        double r = glm::length(r_vect);

        if (p1 != p2) {
            return -spikyGradient(r_vect, r) / (rest_density);
        }

        glm::dvec2 res = glm::dvec2();

        for (int i = 0; i < neighbors[idx1].size(); i++) {
            // Compute difference
            r_vect = p1->ep - particles->at(neighbors[idx1][i])->ep;
            r = glm::length(r_vect);
            res += spikyGradient(r_vect, r);
        }

        return res / (rest_density);
    }

    void updateN(int *sim_n) {}

    double rest_density;
    int num_Particles;
    bool open;

    std::vector<int> particle_indices;
    std::vector<int> *neighbors;

    // This is the same lambda and delta values as explained in Muller's Position-based Fluids
    glm::dvec2 *delta_values;
    std::map<int, double> lambda_values; 
};