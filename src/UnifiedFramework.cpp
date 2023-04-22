#include "UnifiedFramework.hpp"
#include "Constraints/ParticleContactConstraint.hpp"
#include "Constraints/BoundaryConstraint.hpp"

UnifiedFramework::UnifiedFramework() {
    sim_n = NULL;
    init(TWO_FLUID);
    //init(GAS_CLOSED);
    //init(GAS_OPEN);
}

UnifiedFramework::~UnifiedFramework() {
    clear();
}

// TODO: Double check
void UnifiedFramework::clear() {

    // Removing all particles
    for(int i = particles.size()-1; i >= 0; i--) {
        Particle *p = particles[i];
        particles.erase(particles.begin() + i);
        delete(p);
    }

    // Removing all gas Injectors
    for(int i = gasInjectors.size()-1; i >= 0; i--) {
        GasInjector *gi = gasInjectors[i];
        gasInjectors.erase(gasInjectors.begin() + i);
        delete(gi);
    }

    for (int i = 0; i < NUM_CONSTRAINT_GROUPS; i++) {
        if(globalConstraints.find((ConstraintGroup) i) != globalConstraints.end()) {
            std::vector<Constraint *> constrintGroup = globalConstraints[(ConstraintGroup) i];
            for (int j = constrintGroup.size()-1; j >=0; j--) {
                Constraint *c = constrintGroup[j];
                // Finally deleting the actual pointer to constaint
                delete(c);
            }
        }
    }
    globalConstraints.clear();

    // Deleting the N parameters for iterative solver
    if (sim_n) {
        delete[] sim_n;
    }
}

void UnifiedFramework::init(Scene type) {
    // Before each change of scene first clear the previous scenes
    this->clear();
    isGridShown = true;

    // For Gas Particles we use a parameter ALPHA to scale gravity
    gravity = glm::dvec2(0,-9.8);

    switch (type) {
    case GAS_OPEN:
        gasOpen(); break;
    case GAS_CLOSED:
        gasClosed(); break;
    case TWO_FLUID:
        twoFluid(); break;
    }

    sim_n = new int[particles.size()];
}

GasConstraint *UnifiedFramework::createGasConstraint(std::vector<Particle*> *in_particles, double density, Boundary boundary) {
    
    int start_index = particles.size();

    std::vector<int> particle_indices;

    for (int i = 0; i < in_particles->size(); i++) {

        Particle *p = (*in_particles)[i];

        // Making sure that there is no point of infinite mass
        if (p->i_mass == 0.0) {
            exit(1);
        }

        p->ph = GAS;

        particles.push_back(p);

        // Saving the indices of gas particles in the global vector for processing the constraint
        particle_indices.push_back(start_index + i);
    }
    
     GasConstraint *constraint;

    if (boundary == OPEN) {
        constraint = new GasConstraint(density, &particle_indices, true);
    } else if (boundary == CLOSED) {
        constraint = new GasConstraint(density, &particle_indices, false);
    }

    // Adding GasConstraint to the map of global constraints
    if (globalConstraints.find(GENERAL) != globalConstraints.end()) {
        // Map Key exists 
        globalConstraints[GENERAL].push_back(constraint);
    } else {
        //Insert key-value pair for the first time to the map
        std::vector<Constraint*> c = {constraint};
        globalConstraints.insert(std::pair<ConstraintGroup, std::vector<Constraint *>>(GENERAL, c));
    }

    return constraint;
}

FluidConstraint *UnifiedFramework::createFluidConstraint(std::vector<Particle *> *in_particles, double density) {
    int start_index = particles.size();

    std::vector<int> particle_indices;

    for (int i = 0; i < in_particles->size(); i++) {
        Particle *p = in_particles->at(i);

        // Setting particle phase to fluid
        p->ph = FLUID;

        // Double checking that the fluid doesn't have infinite mass
        if (p->i_mass == 0.0) {
            exit(1);
        }
        // updating global particles vector
        particles.push_back(p);

        // updating the local particle_indices for Constraints
        particle_indices.push_back(start_index + i);
    }

    FluidConstraint *constraint = new FluidConstraint(density, &particle_indices);

    // Adding Fluid Constraint to the map of global constraints
    if (globalConstraints.find(GENERAL) != globalConstraints.end()) {
        // Map Key exists 
        globalConstraints[GENERAL].push_back(constraint);
    } else {
        // Insert key-value pair for the first time to the map
        std::vector<Constraint*> c = {constraint};
        globalConstraints.insert(std::pair<ConstraintGroup, std::vector<Constraint *>>(GENERAL, c));
    }

    return constraint;
}

void UnifiedFramework::addGasInjector(glm::dvec2 position, GasConstraint *constraint,  double particlesPerSec, bool isOpen) {
    gasInjectors.push_back(new GasInjector(position, constraint, particlesPerSec, isOpen));
}

void UnifiedFramework::gasClosed() {
    double dx = 0.75;

    xBnd = glm::dvec2(-6.0, 6.0);
    yBnd = glm::dvec2(-6.0, 6.0);

    // Creating Fluid Particles
    std::vector<Particle*> particles;
    glm::vec3 particle_color = glm::vec3(0.36, 0.25f, 0.82f); 

    for(double x = -6.0; x < 6.0; x += dx) {
        for(double y = -6.0; y < 6.0; y += dx) {
            // Add a bit of random extra to the coordinates to make it more realistic
            float random_extra = (((float) rand() / RAND_MAX) * 0.2);
            particles.push_back(new Particle(x+random_extra, y+random_extra, 1, particle_color));
        }
    }

    // Creating a gas constraint
    GasConstraint *constraint = createGasConstraint(&particles, 1.5, CLOSED);

    // Gas Injector explained in section 7.2.1
    addGasInjector(glm::dvec2(0,-2), constraint, 15, false);

    particles.clear();
}

void UnifiedFramework::gasOpen() {
    double dx = 0.6;

    xBnd = glm::dvec2(-8, 8);
    yBnd = glm::dvec2(-6, 200);

    std::vector<Particle*> particles;
    glm::vec3 particle_color = glm::vec3(0.46, 0.85f, 0.82f); 

    for(double x = -6.; x <  6.; x += dx) {
        for(double y = -4.; y < 4.; y += dx) {
            // Add a bit of random extra to the coordinates to make it more realistic
            float random_extra = (((float) rand() / RAND_MAX) * 0.2);
            particles.push_back(new Particle(x+random_extra, y+random_extra, 1, particle_color));
        }
    }

    GasConstraint *constraint = createGasConstraint(&particles, 1.5, OPEN);
    addGasInjector(glm::dvec2(0,-3.), constraint, 20, true);

    // Free mem
    particles.clear();
}

void UnifiedFramework::twoFluid(){
    double dx = 0.62;

    // Setting boundaries
    xBnd = glm::dvec2(-6.0, 6.0);
    yBnd = glm::dvec2(-8.0, 30.0);

    // Creating particles
    std::vector<Particle *> particles;

    double num_fluids = 2.0;
    glm::vec3 particle_color;

    for (int d = 0; d < num_fluids; d++) {
        if (d == 0) particle_color = glm::vec3(0.30, 0.73f, 0.09f);
        else if (d == 1) particle_color = glm::vec3(0.0f, 0.27f, 0.67f);

        double start = -6 + 12 * (d / num_fluids);
        for(double x = start; x < start + (12 / num_fluids); x += dx) {
            for(double y = -8.0; y < 4.0; y += dx) {
                // Add a bit of random extra to the coordinates to make it more realistic
                float random_extra = (((float) rand() / RAND_MAX) * 0.2); 
                particles.push_back(new Particle(x+random_extra, y+random_extra, 1, particle_color, FLUID));
            }
        }

        // Each fluid will be created with different density values
        if (d == 0) {
            createFluidConstraint(&particles, 1);
        } else {
            createFluidConstraint(&particles, 2);
        }

        particles.clear();
    }
}

void UnifiedFramework::showGrid() {
    isGridShown = true;
}

void UnifiedFramework::unshowGrid() {
    isGridShown = false;
}

void UnifiedFramework::display() {
    if (isGridShown)
        displayGrid();
    displayBoundaries();
    displayParticles();
    fixContactConstraints();
    displaySmoke();
}

void UnifiedFramework::resize(int width, int height) {
    // Scale of the GL Ortho is set to 15 by default
    dimensions = glm::ivec2(15,15);
}

void UnifiedFramework::displayGrid() {

    glColor3f(0.3, 0.3, 0.3);
    glBegin(GL_LINES);

    // Grid Vertical Lines 
    for (int x = -dimensions.x; x <= dimensions.x; x++) {
        glVertex2f(x, -dimensions.y);
        glVertex2f(x, dimensions.y);
    }

    // Grid Horizontal Lines
    for (int y = -dimensions.y; y <= dimensions.y; y++) {
        glVertex2f(-dimensions.y, y);
        glVertex2f(dimensions.y, y);
    }

    glColor3f(1,1,1);

    glVertex2f(-dimensions.x, 0);
    glVertex2f(dimensions.x, 0);
    glVertex2f(0, -dimensions.y);
    glVertex2f(0, dimensions.y);

    glEnd();
}

void UnifiedFramework::displayBoundaries() {
    // drawing the boundaries
    glLineWidth(3);

    glBegin(GL_LINES);
    glVertex2f(xBnd.x, yBnd.x);
    glVertex2f(xBnd.x, yBnd.y);

    glVertex2f(xBnd.y, yBnd.x);
    glVertex2f(xBnd.y, yBnd.y);

    glVertex2f(xBnd.x, yBnd.x);
    glVertex2f(xBnd.y, yBnd.x);

    glVertex2f(xBnd.x, yBnd.y);
    glVertex2f(xBnd.y, yBnd.y);
    glEnd();

    glLineWidth(1);
}

void UnifiedFramework::displayParticles() {

    for (int i = 0; i < particles.size(); i++) {
        const Particle *p = particles[i];

        if (p->i_mass == 0.f) {
            glColor3f(1,0,0);
        } else {
            glColor3f(p->color[0], p->color[1], p->color[2]);
        }

        glPushMatrix();
        glTranslatef(p->p.x, p->p.y, 0);
        glScalef(PARTICLE_RAD, PARTICLE_RAD, 0);

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0,0);
        for (int f = 0; f <= 32; f++) {
            double a = f * M_PI / 16.f;
            glVertex2f(sin(a), cos(a));
        }
        glEnd();

        glPopMatrix();
    }
}

void UnifiedFramework::fixContactConstraints() {
    for (int i = 0; i < globalConstraints.size(); i++) {
        if ((ConstraintGroup) i == CONTACT) {
            for (int j = 0; j < globalConstraints[(ConstraintGroup) i].size(); j++) {
                // Fixing particle location that are in contact with a body
                dynamic_cast<ParticleContactConstraint*>(globalConstraints[(ConstraintGroup)i ][j])->fixContactDisplay(&particles);
            }
        }
    }
}

int UnifiedFramework::getNumParticles() {
    return particles.size();
}

// It's the same as fluid particles (circles) but with a smaller radius 
void UnifiedFramework::displaySmoke(){

    glBegin(GL_QUADS);

    double r = PARTICLE_RAD/4;

    for (int i = 0; i < gasInjectors.size(); i++) {
        std::vector<Particle *> *smoke_particles = gasInjectors[i]->getParticles();
        int size = smoke_particles->size();

        for(int j = 0; j < size; j++) {
            Particle *p = (*smoke_particles)[j];
            glColor3f(p->color[0], p->color[1], p->color[2]);
            //cout << "smoke particles at (" << p->p.x << "," << p->p.y << ")" << endl;
            glVertex2d(p->p.x-r, p->p.y-r);
            glVertex2d(p->p.x+r, p->p.y-r);
            glVertex2d(p->p.x+r, p->p.y+r);
            glVertex2d(p->p.x-r, p->p.y+r);
        }
    }
    glEnd();
}

// This is as per Alg 1 of the Muller paper
void UnifiedFramework::advanceTime(double elapsed) {

    // (Step 1) For all Particles apply forces
    for (int i = 0; i < particles.size(); i++) {
        Particle *p = particles[i];

        // For gasses: Consider Alpha * gravity
        // Update velocity and focrce based on f_ext
        if(p->ph == GAS) {
            p->v = p->v + elapsed * (gravity * GRAVITY_REDUCTION_ALPHA) + elapsed * p->f;
        } else {
            p->v = p->v + elapsed * gravity + elapsed * p->f;
        }

        // clearing force vector of p 
        p->clearForce();

        // Predict Position 
        if (p->i_mass == 0) {   // if static
            p->ep = p->p;
        } else {
            p->ep = p->p + p->v * elapsed;
        }

        // Counts for iterative particle solver
        sim_n[i] = 0;

        // Apply mass scaling [m i = mie^(-kh(x_i))] Eq 21 of the paper
        p->scaleMass();
    }

    // Local copy of all constraints to process 
    std::map<ConstraintGroup, std::vector<Constraint *> > constraints;

    for (int i = 0; i < globalConstraints.size(); i++) {
        std::vector<Constraint *> constrGroup = globalConstraints[(ConstraintGroup) i];
        for (int j = 0; j < constrGroup.size(); j++) {
            constraints[(ConstraintGroup) i].push_back(constrGroup[j]);
        }
    }
    
    // (Step 2) For all particles
    for (int i = 0; i < particles.size(); i++) {
        Particle *p = particles[i];

        // Find solid contacts
        for (int j = i + 1; j < particles.size(); j++) {
            Particle *p2 = particles[j];

            // Corner case: Particles that don't move
            if (p->i_mass == 0 && p2->i_mass == 0) {
                continue;
            } else {
                // Detect collision based on particles' circle radius
                if (glm::distance(p->ep, p2->ep) < PARTICLE_DIAMETER - EPSILON) {
                    // Check for Solid-other contact
                    if (p->ph == RIGID || p2->ph == RIGID) {
                        constraints[CONTACT].push_back(new ParticleContactConstraint(i, j));
                    }
                }
            }
        }

        // Boundary Constraints
        if (p->ep.x > xBnd.y - PARTICLE_RAD) {
            constraints[CONTACT].push_back(new BoundaryConstraint(i, xBnd.y, 'x', 'l'));
        } else if (p->ep.x < xBnd.x + PARTICLE_RAD) {
            constraints[CONTACT].push_back(new BoundaryConstraint(i, xBnd.x, 'x', 'g'));
        } 

        if (p->ep.y > yBnd.y - PARTICLE_RAD) {
            constraints[CONTACT].push_back(new BoundaryConstraint(i, yBnd.y, 'y', 'l'));
        } else if (p->ep.y < yBnd.x + PARTICLE_RAD) {
            constraints[CONTACT].push_back(new BoundaryConstraint(i, yBnd.x, 'y', 'g'));
        } 
    }

    // (Step 3) For each constraint group, solve them

    //  Updating n for each constraint group
    for (int j = 0; j < (int) NUM_CONSTRAINT_GROUPS; j++) {
        ConstraintGroup group = (ConstraintGroup) j;

        // Checking the added contact constraints
        //cout << " Constraint is " << group << " and size of it is " << constraints[group].size() << endl;
        for (int k = 0; k < constraints[group].size(); k++) {
            constraints[group].at(k)->updateN(sim_n);
        }
    }

    // Calling project for each constraint 
    for (int i = 0; i < SOLVER_ITRS; i++) {
        for (int j = 0; j < (int) NUM_CONSTRAINT_GROUPS; j++) {
            ConstraintGroup group = (ConstraintGroup) j;
            // Solve constraints and update position
            for (int k = 0; k < constraints[group].size(); k++) {
                constraints[group].at(k)->project(&particles, sim_n);
            }
        }
    }

    // Updating particle velocities [including all the delta_x applied by constrints]
    for (int i = 0; i < particles.size(); i++) {
        Particle *p = particles[i];
        p->v = (p->ep - p->p) / elapsed;
        p->update_pos_epsilon_check();
    }
  
    // Clearing contact constraints for next iterations! 
    for(int i = constraints[CONTACT].size()-1; i >= 0; i--) {
        Constraint *c = constraints[CONTACT][i];
        constraints[CONTACT].erase((constraints[CONTACT]).begin() + i);
        delete(c);
    }

    // Also call advanceTime on injectors to emit the new particles
    for(GasInjector *gi: gasInjectors) {
        gi->advanceTime(&particles, elapsed);

        // Need to check Boundary constraints for the new particles
        for(Particle *p: *(gi->getParticles())) {

            // check xBnd
            if (p->p.x > xBnd.y) {
                p->p.x = xBnd.y;
            } else if (p->p.x < xBnd.x) {
                p->p.x = xBnd.x;
            } 

            // check yBnd
            if (p->p.y > yBnd.y) {
                p->p.y = yBnd.y;
            } else if (p->p.y < yBnd.x) {
                p->p.y = yBnd.x;
            } 
        }
    }

    // clearing parameter n for next iteration
    delete[] sim_n;
    sim_n = new int[particles.size()];
}