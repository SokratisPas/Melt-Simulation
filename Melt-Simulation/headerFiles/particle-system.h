#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>


// =======================================
// ========== Particle Class =============

struct Particle {
    glm::vec3 position = glm::vec3(0.0f);       // possition
    glm::vec3 velocity = glm::vec3(0.0f);       // velocity
    glm::vec3 acceleration = glm::vec3(0.0f);   // acceleration

    glm::vec3 force = glm::vec3(0.0f);          // force

    glm::vec3 accelerationPrev = glm::vec3(0.0f);   // previous acceleration

    float mass = 1.0f;                          // mass (lj)
};


// =======================================
// ========== System Class ===============

class ParticleSystem {
public:
    int partNumber = 1;         // Number of particles

    float Ekin = 0.0f;          // Kinetik energy
    float Epot = 0.0f;          // potential energy
    float temperature = 0.0f;   // temp
    float pressure = 0.0f;      // Pressure
    float virial = 0.0f;        // Virial Pot Energy

    float volume = 1.0f;        // volume
    float density = 1.0f;       // Density

    float dof = 1.0f;           // degrees of Freedom

    // Cube dimensions 
    float xMin = 0.0f, xMax = 1.0f;
    float yMin = 0.0f, yMax = 1.0f;
    float zMin = 0.0f, zMax = 1.0f;

    // vector of particles
    std::vector<Particle> particles;

    // each particle[i] has a vector that includes all the neighbor particles
    std::vector<std::vector<int>> neighborList;

    // last pos of particle[i] 
    std::vector<glm::vec3> lastBuildPos;

    // ===============================================
    // ================== Functions ==================
    // 
    // ----------------------------------------
    // Constructor
    ParticleSystem(int in_partNumber)
    {
        partNumber = in_partNumber;

        // degrees of freedom (we remove 3 bcz we remove net momentum)
        dof = 3.0f * partNumber - 3.0f;

        particles.resize(partNumber);

        neighborList.resize(partNumber);

        lastBuildPos.resize(partNumber);
    }

    // ----------------------------------------
    // calc random float in [start, finish]
    float randomFloat(float start, float finish)
    {
        static std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<float> dist(start, finish);
        return dist(gen);
    }

    // ----------------------------------------
    // particles positions in FCC latice
    std::vector<std::array<double, 3>> generateFCC(int n, double a)
    {
        std::vector<std::array<double, 3>> positions;

        // FCC has 4 atoms per unit cell
        int cells_per_dim = std::ceil(std::cbrt(n / 4.0));

        // Basis positions inside one FCC cell
        std::vector<std::array<double, 3>> basis = {
            {0.0, 0.0, 0.0},
            {0.0, 0.5, 0.5},
            {0.5, 0.0, 0.5},
            {0.5, 0.5, 0.0}
        };

        for (int i = 0; i < cells_per_dim; i++) {
            for (int j = 0; j < cells_per_dim; j++) {
                for (int k = 0; k < cells_per_dim; k++) {

                    for (const auto& b : basis) {
                        if (positions.size() >= n) return positions;

                        positions.push_back({
                            a * (i + b[0]),
                            a * (j + b[1]),
                            a * (k + b[2])
                            });
                    }
                }
            }
        }
        return positions;
    }

    // ---------------------------------------------------
    // initialize Pos, Vel, lastBuildPos (NOT Acceleration)
    void initializeParameteres(double laticeConst) 
    {
        std::vector<std::array<double, 3>> initialPos = generateFCC(partNumber, laticeConst); // initial pos vector

        for (int i = 0; i < partNumber; i++)
        {
            // initial pos
            particles[i].position = glm::vec3(initialPos[i][0], initialPos[i][1], initialPos[i][2]);

            // initial small random vel
            particles[i].velocity = glm::vec3(randomFloat(-0.1f, 0.1f), randomFloat(-0.1f, 0.1f),
                randomFloat(-0.1f, 0.1f));

            // initial lastBuildPos
            lastBuildPos[i] = particles[i].position;
        }
    }

    // ----------------------------------------
    // remove net momentum
    // it should : sum(vel[i] == 0)
    void removeNetMomentum()
    {
        glm::vec3 velMean(0.0f);
        for (auto& p : particles)
            velMean += p.velocity;

        velMean /= (float)particles.size();

        for (auto& p : particles)
            p.velocity -= velMean;
    }
      
    // ----------------------------------------
    // aplly PBC 
    void applyPBC(glm::vec3& pos)
    {
        // cube lengths
        float Lx = xMax - xMin;
        float Ly = yMax - yMin;
        float Lz = zMax - zMin;

        // keep x, y, x in cube boundaries 
        pos.x = xMin + fmod(fmod(pos.x - xMin, Lx) + Lx, Lx);
        pos.y = yMin + fmod(fmod(pos.y - yMin, Ly) + Ly, Ly);
        pos.z = zMin + fmod(fmod(pos.z - zMin, Lz) + Lz, Lz);
    }

    // ----------------------------------------
    // Build Neighbour List
    void buildNeighborList(float rlist)
    {
        float rlist2 = rlist * rlist;

        int N = partNumber;

        // clear neighborList
        for (auto& nl : neighborList)
            nl.clear();

        // build cells (using rlist as cell size)
        float cellSize = rlist;

        // number of cels
        int nx = (int)((xMax - xMin) / cellSize);
        int ny = nx;
        int nz = nx;

        // debug line for n <= 0
        if (nx <= 0 || ny <= 0 || nz <= 0)
        {
            std::cout << "Invalid cell grid!\n";
            return;
        }

        std::vector<std::vector<int>> cells(nx * ny * nz);

        // make 1D coords to 3D coords
        auto getCellIndex = [&](int ix, int iy, int iz)
            {
                return ix + nx * (iy + ny * iz);
            };

        // assign particles to cells
        for (int i = 0; i < N; i++)
        {
            glm::vec3 p = particles[i].position;

            // cells indices
            int ix = (int)((p.x - xMin) / cellSize) % nx;
            int iy = (int)((p.y - yMin) / cellSize) % ny;
            int iz = (int)((p.z - zMin) / cellSize) % nz;

            // catch negative indices
            if (ix < 0) ix += nx;
            if (iy < 0) iy += ny;
            if (iz < 0) iz += nz;

            cells[getCellIndex(ix, iy, iz)].push_back(i);
        }

        // build neighbor list
        for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
                for (int iz = 0; iz < nz; iz++)
                {
                    int cellA = getCellIndex(ix, iy, iz);

                    // check only neighbors of cellA
                    for (int dx = -1; dx <= 1; dx++)
                        for (int dy = -1; dy <= 1; dy++)
                            for (int dz = -1; dz <= 1; dz++)
                            {
                                int jx = (ix + dx + nx) % nx;
                                int jy = (iy + dy + ny) % ny;
                                int jz = (iz + dz + nz) % nz;

                                int cellB = getCellIndex(jx, jy, jz);

                                for (int i : cells[cellA])
                                {
                                    for (int j : cells[cellB])
                                    {
                                        // avoid double counting
                                        if (i >= j) continue;

                                        glm::vec3 r = particles[i].position - particles[j].position;

                                        // minimum image
                                        float Lx = xMax - xMin;
                                        float Ly = yMax - yMin;
                                        float Lz = zMax - zMin;

                                        r.x -= Lx * std::round(r.x / Lx);
                                        r.y -= Ly * std::round(r.y / Ly);
                                        r.z -= Lz * std::round(r.z / Lz);

                                        float dist2 = glm::dot(r, r);

                                        if (dist2 < rlist2)
                                        {
                                            neighborList[i].push_back(j);
                                            neighborList[j].push_back(i);
                                        }
                                    }
                                }
                            }
                }

        // store positions
        for (int i = 0; i < N; i++)
            lastBuildPos[i] = particles[i].position;
    }

    // ----------------------------------------
    // particle-particle lj force
    glm::vec3 computeLJForce(const glm::vec3& pi, const glm::vec3& pj, float epsilon, float sigma)
    {
        glm::vec3 r = pi - pj; // distance vec of pair

        // box dimensions
        float Lx = xMax - xMin;
        float Ly = yMax - yMin;
        float Lz = zMax - zMin;

        // minimum image convention
        // IDEA:
        // use the smallest distance between 2 particles 
        // r / L callc how many boxes is the distance 
        // and we remove this distance from r
        r.x -= Lx * std::round(r.x / Lx);
        r.y -= Ly * std::round(r.y / Ly);
        r.z -= Lz * std::round(r.z / Lz);

        float dist2 = glm::dot(r, r);

        // avoid division by zero 
        if (dist2 < 1e-12f)
            return glm::vec3(0.0f);

        // cut off radius
        float rc = 2.5f * sigma;
        if (dist2 > (rc * rc))
            return glm::vec3(0.0f);

        // distance magnitude
        float dist = std::sqrt(dist2);

        // inverse distance
        float inv_r = glm::inversesqrt(dist2);

        float s_over_r = sigma * inv_r;
        float s2 = s_over_r * s_over_r;
        float s6 = s2 * s2 * s2;
        float s12 = s6 * s6;

        // E potential update
        Epot += 4.0f * epsilon * (s12 - s6);

        // LJ force magnitude
        float f = 24.0f * epsilon * (2.0f * s12 - s6) * inv_r;

        return r * (f * inv_r);
    }

    // ----------------------------------------
    // calc all lj forces and Epot and Virial 
    void computeForces_Epot_Virial(float epsilon, float sigma)
    {
        int N = partNumber;

        // set forces = 0
        for (auto& p : particles)
            p.force = glm::vec3(0.0f);

        // calc new lj forces
        for (int i = 0; i < N; i++)
        {
            for (int j : neighborList[i])
            {
                if (i >= j) continue;

                glm::vec3 ijForce = computeLJForce(particles[i].position, particles[j].position, epsilon, sigma);

                particles[i].force += ijForce;
                particles[j].force -= ijForce;


                glm::vec3 rij = particles[i].position - particles[j].position;
                virial += glm::dot(rij, ijForce);
            }
        }
    }

    // ----------------------------------------
    // calc volume, density
    void calcVolumeDensity()
    {
        float volume = (xMax - xMin) * (yMax - yMin) * (zMax - zMin);
        float density = partNumber / volume;
    }

    // ----------------------------------------
    // initiallize acceleration
    void initiallizeAcc()
    {
        for (auto& p : particles)
        {
            p.acceleration = p.force;       // m = 1
            p.accelerationPrev = p.force;   // m = 1
        }
    }

    // ----------------------------------------
    // reset Ekin, Epot, Virial
    void resetEkinEpotVirial()
    {
        Ekin = 0.0f, Epot = 0.0f, virial = 0.0f;
    }

    // ----------------------------------------
    // update pos wth velocity Verlet method
    void updatePos(float dt)
    {
        for (auto& part : particles)
        {
            part.position += part.velocity * (float)dt +
                0.5f * part.accelerationPrev * (float)(dt * dt);

            // periodic boundary conditions
            applyPBC(part.position);

            // save prev acceleration
            part.accelerationPrev = part.acceleration;
        }
    }

    // ----------------------------------------
    // check if neighboor list needs rebuild
    // IDEA: If any particle has moved (distance > rs / 2), we rebuild the neighborList
    bool needRebuild(float skin)
    {
        int N = partNumber;

        // max displacement^2 of any particle
        float maxDisp2 = 0.0f;

        for (int i = 0; i < N; i++)
        {
            glm::vec3 dr = particles[i].position - lastBuildPos[i];

            // minimum image
            float Lx = xMax - xMin;
            float Ly = yMax - yMin;
            float Lz = zMax - zMin;

            dr.x -= Lx * std::round(dr.x / Lx);
            dr.y -= Ly * std::round(dr.y / Ly);
            dr.z -= Lz * std::round(dr.z / Lz);

            maxDisp2 = std::max(maxDisp2, glm::dot(dr, dr));
        }

        return maxDisp2 > (skin * 0.5f) * (skin * 0.5f);
    }

    // ----------------------------------------
    // acceleration update
    void updateAcc()
    {
        for (auto& part : particles)
        {
            part.acceleration = part.force; // m = 1
        }
    }

    // ----------------------------------------
    // vel and Ekin update
    void updateVelEkin(float dt)
    {
        for (auto& part : particles)
        {
            part.velocity += 0.5f * (part.accelerationPrev + part.acceleration) * (float)dt;

            Ekin += 0.5f * glm::dot(part.velocity, part.velocity);
        }
    }

    // ----------------------------------------
    // calc T (kb = 1, lj units)
    void computeTemperature()
    {
        temperature = (2.0f * Ekin) / dof;
    }

    // ----------------------------------------
    // calc pressure (= density * temp + sum(Rij*Fij) / (3V) )
    void computePressure()
    {
        pressure =  density * temperature + virial / (3.0f * volume);
    }

    // ----------------------------------------
    // thermostat (Berendsen)
    void useThermostat(const float dt, const float tau, const float Ttarget)
    {
        if (temperature < 1e-8f) return; // avoid division by zero

        // scaling factor
        float lambda = sqrt(1.0f + (dt / tau) * (Ttarget / temperature - 1.0f));

        // change velocity of particles
        for (auto& p : particles)
            p.velocity *= lambda;
    }
};