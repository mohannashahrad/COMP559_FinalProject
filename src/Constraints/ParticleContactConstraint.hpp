#include "../Particle.hpp"

// Constraint for contact between two particles that at least one of them is a rigid body
class ParticleContactConstraint : public Constraint
{
public:
    ParticleContactConstraint(int p1, int p2, bool is_static = false) : Constraint() {
        particle1 = p1;
        particle2 = p2;
        isStatic = is_static;
    }

    virtual ~ParticleContactConstraint(){}

    // Particles have the estimated value
    void project(std::vector<Particle *> *particles, int *sim_n) {

        Particle *p1 = (*particles)[particle1];
        Particle *p2 = (*particles)[particle2];

        // The scaled mass (Eq 21 of paper) is used for contact constraints
        if (p1->t_mass == 0.f && p2->t_mass == 0.f) {
            return;
        }


        glm::dvec2 diff;
        if (isStatic) {
            diff = p1->p - p2->p;
        } else {
            diff = p1->ep - p2->ep;
        }

        double distance = glm::length(diff);

        // There is no collision
        if (distance > PARTICLE_DIAMETER) {
            return;
        }

        // Case where there is a collision between the particles

        glm::dvec2 delta_x = diff * ((distance - PARTICLE_DIAMETER) / (p1->t_mass + p2->t_mass) / distance);
        glm::dvec2 delta_x1 = -p1->t_mass * delta_x / (double)sim_n[particle1];
        glm::dvec2 delta_x2 = p2->t_mass * delta_x / (double)sim_n[particle2];

        // Updating position based on the displacements
        p1->ep += delta_x1;
        p2->ep += delta_x2;

        // If particles are not moving
        if (isStatic) {
            p1->p += delta_x1;
            p2->p += delta_x2;
        }
    }

    void fixContactDisplay(std::vector<Particle *> *particles) {

        Particle *p1 = (*particles)[particle1];
        Particle *p2 = (*particles)[particle2];

        glColor3f(1,1,0);
        glBegin(GL_LINES);

        glVertex2f(p1->p.x, p1->p.y);
        glVertex2f(p2->p.x, p2->p.y);

        glEnd();

        glPointSize(3);
        glBegin(GL_POINTS);

        glVertex2f(p1->p.x, p1->p.y);
        glVertex2f(p2->p.x, p2->p.y);

        glEnd();
    }
    
    void updateN(int *sim_n) {
        sim_n[particle1]++;
        sim_n[particle2]++;
    }

private:
    int particle1, particle2;
    bool isStatic;
};
