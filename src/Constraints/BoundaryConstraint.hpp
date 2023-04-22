#include "../Particle.hpp"

// This constraint is for the collisions with boundaries
class BoundaryConstraint : public Constraint
{
public:
    BoundaryConstraint(int i, double v, char bndType, char cmp, bool is_static = false) : Constraint() {
        index  = i;
        limit = v;
        boundaryType = bndType;
        comparison = cmp;
        isStatic = is_static;
    }
    virtual ~BoundaryConstraint() {}

    // Look at section 6.1 of the unified particle system paper for this
    void project(std::vector<Particle *> *particles, int *sim_n) {
        
        // Nothing to do if the particle is isStatic
        if (isStatic) {
            return;
        }

        glm::dvec2 normal = glm::dvec2();
        Particle *p = particles->at(index);

        // For fluid and gasses around boundaries
        double randomlimit = ((double)rand() / (double)RAND_MAX) * 0.0025;

        double extra_rad = p->ph == FLUID || p->ph == GAS ? randomlimit : 0;
        double d = (PARTICLE_RAD + extra_rad);

        // Move particles back into a valid spot (if passes the boundaries)
        if (comparison == 'g' && boundaryType == 'x') {
            if (p->ep.x >= limit + PARTICLE_RAD) return;

            normal = glm::dvec2(1,0);
            p->ep.x = limit + d;

            if (isStatic) 
                p->p.x = limit + d;

        } else if (comparison == 'g' && boundaryType == 'y') {
            if (p->ep.y >= limit + PARTICLE_RAD) return;

            normal = glm::dvec2(0,1);
            p->ep.y = limit + d;
            if (isStatic) 
                p->p.y = limit + d;

        } else if (comparison == 'l' && boundaryType == 'x') {
            if (p->ep.x <= limit - PARTICLE_RAD) return;

            normal = glm::dvec2(-1,0);
            p->ep.x = limit - d;
            if (isStatic) 
                p->p.x = limit - d;
        
        } else if (comparison == 'l' && boundaryType == 'y') {
            if (p->ep.y <= limit - PARTICLE_RAD) return;

            normal = glm::dvec2(0,-1);
            p->ep.y = limit - d;
            if (isStatic) 
                p->p.y = limit - d;
        }

        // Boundaries have a coefficient of friction of 1
        glm::dvec2 dp = (p->ep - p->p) / (double)sim_n[index];
        glm::dvec2 dpt = dp - glm::dot(dp, normal) * normal;

        double l = glm::length(dpt);

        // Epsilon check 
        if (l < EPSILON) {
            return;
        }

        // Static and kinetic friction
        // Equation(24) of the Muller Unified Particle System Paper
        if (l < sqrt(p->sFriction) * d) {
            p->ep -= dpt;
        } else {
            p->ep -= dpt * min(sqrt(p->kFriction )* d / l, 1.0);
        }
    }

    void updateN(int *sim_n) {
        sim_n[index]++;
    }

private:
    int index;
    double limit;
    bool isStatic;
    char boundaryType;
    char comparison;
};
