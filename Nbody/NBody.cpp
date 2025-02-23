// #include <iostream>
// #include <vector>
// #include <string>

// using namespace std;

// struct Particle{
//     double mass;
//     double position[3];

// }

// int main(){
//     int particleNo = 3;

//     double masses[particleNo];
//     double positions[particleNo];
//     double velocities[particleNo];
//     double forces[particleNo];

//     for(int i = 0; i < particleNo; i++){
//             masses[i] = rand();
//             positions[i] = rand();
//             velocities[i] = rand();
//     }

//     return 0;
// }
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

const double G = 6.674e-11;   // Gravitational constant
const double SOFTENING = 0.0001; // Small value to prevent division by zero

// Structure to represent a particle (like a planet or star)
struct Particle {
    double mass;  // Mass of the particle
    double x, y, z; // Position in 3D space
    double vx, vy, vz; // Velocity in 3D space
    double fx, fy, fz; // Force acting on the particle
};

// Function to calculate forces between all pairs of particles
void computeForces(vector<Particle>& particles) {
    int n = particles.size();
    
    // Reset forces to zero before calculating new ones
    for (int i = 0; i < n; i++) {
        particles[i].fx = particles[i].fy = particles[i].fz = 0.0;

        for (int j = 0; j < n; j++) {
            if (i != j) {  // Don't calculate force on itself
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;

                double r2 = dx * dx + dy * dy + dz * dz + SOFTENING; // Distance squared
                double r = sqrt(r2); // Actual distance
                double F = (G * particles[i].mass * particles[j].mass) / r2; // Newton’s Law

                // Add forces in the correct directions
                particles[i].fx += F * dx / r;
                particles[i].fy += F * dy / r;
                particles[i].fz += F * dz / r;
            }
        }
    }
}

// Function to update the positions and velocities of all particles
void updateParticles(vector<Particle>& particles, double dt) {
    for (auto& p : particles) {
        double ax = p.fx / p.mass; // Acceleration in x
        double ay = p.fy / p.mass; // Acceleration in y
        double az = p.fz / p.mass; // Acceleration in z

        // Update velocity using acceleration
        p.vx += ax * dt;
        p.vy += ay * dt;
        p.vz += az * dt;

        // Update position using velocity
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }
}

// Function to save the state of particles to a file
void saveState(const vector<Particle>& particles, ofstream& file) {
    file << particles.size(); // Write number of particles
    for (const auto& p : particles) {
        file << "\t" << p.mass << "\t" << p.x << "\t" << p.y << "\t" << p.z
             << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz
             << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
    }
    file << "\n"; // New line for each step
}

int main(int argc, char* argv[]){
    // Check if the user provided enough arguments
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <num_particles> <dt> <steps> <output_frequency>" << endl;
        return 1;
    }
    
    // Read input values from command line
    int numParticles = atoi(argv[1]); // Number of particles
    double dt = atof(argv[2]); // Time step
    int steps = atoi(argv[3]); // Number of simulation steps
    int outputFreq = atoi(argv[4]); // How often to save data

    // Create particles with random positions and masses
    vector<Particle> particles(numParticles);
    for (auto& p : particles) {
        p.mass = rand() % 100 + 1; // Random mass between 1 and 100
        p.x = rand() % 200 - 100; // Random position between -100 and 100
        p.y = rand() % 200 - 100;
        p.z = rand() % 200 - 100;
        p.vx = p.vy = p.vz = 0; // Start with no velocity
    }

    ofstream file("output.tsv"); // Open file for saving results

    // Run the simulation for the given number of steps
    for (int step = 0; step < steps; step++) {
        computeForces(particles); // Calculate forces
        updateParticles(particles, dt); // Move particles
        
        // Save data every 'outputFreq' steps
        if (step % outputFreq == 0) {
            saveState(particles, file);
        }
    }

    file.close(); // Close the file
    cout << "Simulation complete. Output saved to output.tsv" << endl;
    return 0;
}
