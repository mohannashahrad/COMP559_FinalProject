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

        // Compute Lambda Values (Equation 11 of Mullers's Position-based Fluids)
        lambda_values.clear();
        for (int k = 0; k < particle_indices.size(); k++) {
            neighbors[k].clear();
            int i = particle_indices[k];

            Particle *p_i = particles->at(i);
            double pi = 0., denom = 0.;

            // Find neighbors
            for (int j = 0; j < particles->size(); j++) {

                // Check if the next particle is actually this particle
                if (j != i) {
                    Particle *p_j = particles->at(j);

                    // Nothing to do for fixed particles
                    if (p_j->i_mass == 0) continue;

                    glm::dvec2 r = p_i->ep - p_j->ep;

                    double rlen2 = glm::dot(r, r);

                    if (rlen2 < H2) {
                        // Found a neighbor! Remember it and add to pi and the gamma denominator
                        neighbors[k].push_back(j);
                        double incr = poly6(rlen2) / p_j->i_mass;
                        if (p_j->ph == RIGID) {
                            incr *= S_SOLID_COUPLING;
                        }
                        pi += incr;

                        glm::dvec2 gr = calculateGrad(particles, k, j);
                        denom += glm::dot(gr, gr);
                    }

                // If it is, cut to the chase
                } else {
                    neighbors[k].push_back(j);
                    pi += poly6(0) / p_i->i_mass;
                }
            }

            glm::dvec2 gr = calculateGrad(particles, k, i);
            denom += glm::dot(gr, gr);

            double p_rat = (pi/rest_density);

            // Enforcing of f_drag [Equation 29 of Paper]
            if(open) {
                //p_i->f += p_i->v * (1.-p_rat) * -50.0;
                p_i->addForce(p_i->v * (1.0-p_rat) * -2.0);
            }

            double relaxation = 0.05;   // For gamma correction
            double lambda = -(p_rat - 1.) / (denom + relaxation);
            lambda_values[i] = lambda;
        }

        // Compute actual delta_values
        for (int k = 0; k < particle_indices.size(); k++) {
            glm::dvec2 delta = glm::dvec2();
            glm::dvec2 f_vort = glm::dvec2();
            int i = particle_indices[k];
            Particle *p_i = particles->at(i);

            for (int x = 0; x < neighbors[k].size(); x++) {
                int j = neighbors[k][x];
                if (i == j) continue;

                Particle *p_j = particles->at(j);
                glm::dvec2 r = p_i->ep - p_j->ep;

                double rlen = glm::length(r);
                glm::dvec2 sg = spikyGradient(r, rlen);

                double lambdaCorr = -K_TERM * pow((poly6(rlen * rlen) / poly6(DELTA_Q_TERM_GAS * DELTA_Q_TERM_GAS * H * H)), N_TERM);
                delta += (lambda_values[i] + lambda_values[j] + lambdaCorr) * sg;


                //  vorticity
                glm::dvec2 gradient = spikyGradient(r, glm::dot(r,r));
                glm::dvec2 w = gradient * p_j->v;
                glm::dvec3 cross = glm::cross(glm::dvec3(0,0,glm::length(w)), glm::dvec3(r.x, r.y, 0));
                f_vort += glm::dvec2(cross.x, cross.y) * poly6(glm::dot(r,r));
            }
            delta_values[k] = (delta / rest_density);
            p_i->f += f_vort;
        }

        for (int k = 0; k < particle_indices.size(); k++) {
            int i = particle_indices[k];
            Particle *p_i = particles->at(i);
            p_i->ep += delta_values[k] / ((double) neighbors[k].size() + sim_n[i]);
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