// meltSimulation.cpp
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"

#include "headerFiles/Shader.h"

// ==========================================================================================================
// ==========================================================================================================
// ====================================== meltSimulation.cpp ================================================
// 
// Simulation of a solid with FCC lattice with Lennard-Jones potential. 
// In order to have full cube with atoms in every position available use: 
// Number of atoms = 4 * n^3 , where n = 1, 2, 3, ...
// For example : N = 256, 500, 864, 1372, 2048, 2916, 4000.
#define N 2916
//
// For the use of the thermostat set USE_THERMOSTAT == 1.
#define USE_THERMOSTAT 1
//
// Keep in mind that the code creates a file to store the results.
// For a test run write the results on a fie "dummy.txt" for example.
// You can define how many steps to truck from the "stepMax" variable.
// File name :
#define FILE_NAME "dummy.txt"
//
// CONTROLS :
// - Use RIGHT CLICK to move the camera. 
// - Start to start the simulations (md calculations).
// - Stop to stop the calculations (also writing on the file stops).
// - Close to end the programm. 
// - ESCAPE also closes the window and ends the programm.
// - Show Cube checkbox to render a cube at the bondaries of the simulation. 
// ==========================================================================================================
// ==========================================================================================================

// ==== particle vertex struct =====
struct particleVertex {
    float x, y, z;      // pos
    float nx, ny, nz;   // normals
};

// ===== particles properties =====
float particleRadius = 0.3f; // scale particles
float mass = 1.0f; // lj units

struct Particle {
    glm::vec3 position = glm::vec3(0.0f);
    glm::vec3 velocity = glm::vec3(0.0f);
    glm::vec3 acceleration = glm::vec3(0.0f);

    glm::vec3 force = glm::vec3(0.0f);

    glm::vec3 accelerationPrev = glm::vec3(0.0f);
};
std::vector<Particle> particles(N); 


// ==== Cube dimensions =====
float xMin = 0.0f, xMax = 1.0f;
float yMin = 0.0f, yMax = 1.0f;
float zMin = 0.0f, zMax = 1.0f;

// ==== Instance data struct =====
struct InstanceData {
    glm::vec3 instancePos;
    glm::vec3 instanceColor;
};
std::vector<InstanceData> instances(N);

// === neighbor list ====
// each particle[i] has a vector that includes all the neighbor particles
std::vector<std::vector<int>> neighborList(N);

// last pos of particle[i] 
std::vector<glm::vec3> lastBuildPos(N);


// Starting screen settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 750;

// ==== orbit camera =====
bool rightMousePressed = false;
float yaw = -90.0f;
float pitch = 0.0f;

float fov = 45.0f; // field of view 

double lastMouseX = 0.0f;
double lastMouseY = 0.0f;

float distanceToTarget = 25.0f;
glm::vec3 target(0.5f, 0.5f, 0.5f); // box center

// time
float deltaTime = 0.0f;	// Time between current frame and last frame
float lastFrame = 0.0f; // Time of last frame

// t, step for file
int step = 0;
float simTime = 0.0;
int stepMax = 2000;

// ===== functions =====
void processInput(GLFWwindow* window);   
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
void generateSphere(float radius, unsigned int sectors, unsigned int stacks, 
                    std::vector<particleVertex>& vertices, std::vector<unsigned int>& indices);
float randomFloat(float start, float finish);
std::vector<std::array<double, 3>> generateFCC(int n, double a);
glm::vec3 computeLJForce(const glm::vec3& pi, const glm::vec3& pj, float epsilon, float sigma, float& Epot);
void buildNeighborList(float rlist);
void computeForces_Epot_Virial(float epsilon, float sigma, float& Epot, float& Virial);
void applyPBC(glm::vec3& pos);
bool needRebuild(float skin);
float computeTemperature(const float Ekin, const float dof);
void useThermostat(const float dt, const float tau, const float Ttarget,
                   const float dof, const float Ekin, const float T);
float computePressure(float density, float temp, float virial, float vol);


int main(void)
{
    // ==== GLFW window initialization ========
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE); // full screen
    glfwWindowHint(GLFW_SAMPLES, 4); // multisample buffer for antialiasing

    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Melt simulation", NULL, NULL);
    if (!window)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    // ======== ImGui initialization =======
    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(window, true);
    ImGui::StyleColorsDark();

    bool showWindow = true;
    bool startSimulation = false;
    bool showCube = false;
    bool showCloseWindow = false;
    float windowTemp = 0.0f, windowEnergy = 0.0f, windowPressure;
    int windowStep = 0;

    // ===== Mouse =======
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);


    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // new width, height for full screen
    int fullWidth, fullHeight;
    glfwGetFramebufferSize(window, &fullWidth, &fullHeight);

    // glew Initialization
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK)
    {
        std::cerr << "Error initializing GLEW: " << glewGetErrorString(err) << std::endl;
        glfwTerminate();
        return -1;
    }

    // =========== SHADERS ==============
    Shader particleShader("C:\\Users\\spast\\source\\repos\\Melt-Simulation\\Melt-Simulation\\vertexShaders\\meltSimulation.vert",
                          "C:\\Users\\spast\\source\\repos\\Melt-Simulation\\Melt-Simulation\\fragmentShaders\\meltSimulation.frag");

    Shader cubeShader("C:\\Users\\spast\\source\\repos\\Melt-Simulation\\Melt-Simulation\\vertexShaders\\meltSimulation2.vert",
                      "C:\\Users\\spast\\source\\repos\\Melt-Simulation\\Melt-Simulation\\fragmentShaders\\meltSimulation2.frag");

    // =================================
    // =========== MODEL DATA ==========
    //==================================
    // =========== Particle =============
    std::vector<particleVertex> partVerts;
    std::vector<unsigned int> partIndices;
    unsigned int partSectors = 18;
    unsigned int partStacks = 12;
    generateSphere(particleRadius, partSectors, partStacks, partVerts, partIndices);

    unsigned int partVBO, partVAO, partEBO;
    glGenVertexArrays(1, &partVAO);
    glBindVertexArray(partVAO);

    glGenBuffers(1, &partVBO);
    glGenBuffers(1, &partEBO);

    glBindBuffer(GL_ARRAY_BUFFER, partVBO);
    glBufferData(GL_ARRAY_BUFFER, partVerts.size() * sizeof(particleVertex), partVerts.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, partEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, partIndices.size() * sizeof(unsigned int), partIndices.data(), GL_STATIC_DRAW);

    // vertex pos
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(particleVertex), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normals 
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(particleVertex), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);


    // =========== cube ============
    float cubeVertices[] = {
        0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f,
        1.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f,
        0.0f, 1.0f, 1.0f
    };
    unsigned int cubeIndices[] = {
        0,1, 1,2, 2,3, 3,0,   // back square
        4,5, 5,6, 6,7, 7,4,   // front square
        0,4, 1,5, 2,6, 3,7    // connecting edges
    };

    unsigned int cubeVAO, cubeVBO, cubeEBO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cubeIndices), cubeIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    // ========= Instance Buffer ==========
    unsigned int instanceVBO;
    glGenBuffers(1, &instanceVBO);

    glBindVertexArray(partVAO);

    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(InstanceData), nullptr, GL_DYNAMIC_DRAW);

    // instancePos
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)0);
    glVertexAttribDivisor(2, 1);

    // instanceColor
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, instanceColor));
    glVertexAttribDivisor(3, 1);
    glBindVertexArray(0);


    // ======================================
    // =========== PHYSICS DATA =============
    double e = 1.0;     // epsilon (lj units)
    double sigma = 1.0; // sigma (lj units)

    double r0 = pow(2.0, 1.0 / 6.0) * sigma;
    double a = sqrt(2.0) * r0;   // lattice constant ( = sqrt(2) * 1.122)

    double dt = 0.0005; 
    float Ekin = 0.0f; // Kinetik energy
    float Epot = 0.0f; // potential energy
    float temperature = 0.0f; // temp
    float pressure = 0.0f;
    float virial = 0.0f; 

    float dof = 3.0f * particles.size() - 3.0f; // degrees of freedom (we remove 3 bcz we remove net momentum)

    // thermostat vars
    float Ttarget = 1.0f;   // target temperature (lj)
    float tau = 0.1f;       // relaxation strength 


    // ==== initialize parameteres ====
    std::vector<std::array<double, 3>> initialPos = generateFCC(N, a); // initial pos vector

    for (int i = 0; i < N; i++)
    {
        // initial pos
        particles[i].position = glm::vec3(initialPos[i][0], initialPos[i][1], initialPos[i][2]);

        // initial small random vel
        particles[i].velocity = glm::vec3(randomFloat(-0.1f, 0.1f), randomFloat(-0.1f, 0.1f),
                                          randomFloat(-0.1f, 0.1f));

        // initial lastBuildPos
        lastBuildPos[i] = particles[i].position;
    }

    // remove net momentum
    // it should : sum(vel[i] == 0)
    glm::vec3 velMean(0.0f);
    for (auto& p : particles)
        velMean += p.velocity;

    velMean /= (float)particles.size();

    for (auto& p : particles)
        p.velocity -= velMean;


    // ==== rescale cube =======
    // (rescale cube BEFORE initialize neighborList and computeForces)
    int cells_per_dim = std::ceil(std::cbrt(N / 4.0));
    float L = cells_per_dim * a; // cubes length

    xMin = yMin = zMin = 0.0f;
    xMax = yMax = zMax = L; 

    // volume, density calculation 
    float volume = (xMax - xMin) * (yMax - yMin) * (zMax - zMin);
    float density = N / volume;


    // ==== neighbor list =====
    float rc = 2.5f * sigma;    // r cutoff = 2.5 * 1.0 = 2.5
    float rSkin = 0.3f;         // r skin
    float rlist = rc + rSkin;   // neighbor list radius

    // apply pbc
    for (auto& p : particles)
        applyPBC(p.position);

    // initalize neighbour list and forces
    buildNeighborList(rlist);
    computeForces_Epot_Virial(e, sigma, Epot, virial);

    // initiallize acceleration
    for (auto& p : particles)
    {
        p.acceleration = p.force; // m = 1
        p.accelerationPrev = p.force;
    }


    // change target for orbital camera
    target = glm::vec3(xMax / 2);
    
    // enable depth test
    glEnable(GL_DEPTH_TEST); 

    // face culling
    glEnable(GL_CULL_FACE);

    // anti aliasing
    glEnable(GL_MULTISAMPLE);
     

    // ======= open file ===========    
    std::ofstream meltFile(FILE_NAME);

    if (!meltFile)
    {
        std::cout << "Failed to open" << FILE_NAME << "file\n";
        return -1;
    }

    meltFile << "step time Ekin Epot Etot Temperature Pressure frameDiff\n";


    // ==============================
    // ======== WHILE LOOP ==========
    while (!glfwWindowShouldClose(window))
    {
        processInput(window);

        ImGui_ImplGlfwGL3_NewFrame();

        glClearColor(0.1f, 0.04f, 0.25f, 1.0f); // background color  
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // calc delaTime (needed for moving camera)
        float currentFrame = static_cast<float>(glfwGetTime()); // time t;
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;


        particleShader.use();
        // projection for particleShader  
        glm::mat4 projection = glm::perspective( glm::radians(fov), (float)fullWidth / (float)fullHeight, 0.1f, 100.0f);
        
        // view for particleShader
        glm::vec3 direction;
        direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        direction.y = sin(glm::radians(pitch));
        direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));

        glm::vec3 cameraPos = target - direction * distanceToTarget;

        glm::mat4 view = glm::lookAt(cameraPos, target, glm::vec3(0.0f, 1.0f, 0.0f));
        particleShader.setMat4("projection", projection);
        particleShader.setMat4("view", view);


        // ======== PHYSICS PARAMETERS UPDATE =========  
        if (startSimulation)
        {
            // reset energy
            Ekin = 0.0f, Epot = 0.0f, virial = 0.0f;

            // pos update
            for (auto& part : particles)
            {
                part.position += part.velocity * (float)dt + 
                                0.5f * part.accelerationPrev * (float)(dt * dt);

                // periodic boundary conditions
                applyPBC(part.position);
                
                // save prev acceleration
                part.accelerationPrev = part.acceleration;
            }

            // check if rebuild neighboorList
            if (needRebuild(rSkin))
            {
                buildNeighborList(rlist);
            }

            // update forces and E potential and Virial
            computeForces_Epot_Virial(e, sigma, Epot, virial);

            // acceleration update
            for (auto& part : particles)
            {
                part.acceleration = part.force; // m = 1
            }

            // vel and Ekin update
            for (auto& part : particles)
            {
                part.velocity += 0.5f * (part.accelerationPrev + part.acceleration) * (float)dt;

                Ekin += 0.5f * glm::dot(part.velocity, part.velocity);
            }

            // Temperature
            temperature = computeTemperature(Ekin, dof);

            // Pressure
            pressure = computePressure(density, temperature, virial, volume);

            // thermostat
            if (USE_THERMOSTAT)
            {
                useThermostat(dt, tau, Ttarget, dof, Ekin, temperature);
            }

            // ======= WRITE FILE ========            
            if (step % 10 == 0)
            {
                float Etot = Ekin + Epot;

                meltFile << step << " "
                    << simTime << " "
                    << Ekin << " "
                    << Epot << " "
                    << Etot << " "
                    << temperature << " "
                    << pressure << " "
                    << deltaTime << "\n";
            }
            step++;
            simTime += dt;

            if (step == stepMax + 1)
                showCloseWindow = true;
        }


        // update instance pos and color
        instances.clear();
        for (auto& part : particles)
        {
            InstanceData data;
            // position
            data.instancePos = part.position;

            // velocity^2
            float u2 = glm::dot(part.velocity, part.velocity);

            // average v^2 from temperature
            float u2Mean = 3.0f * temperature;

            // choose upper visualization limit
            float u2Max = 4.0f * u2Mean;

            // avoid divide by zero
            if (u2Max < 1e-6f)
                u2Max = 1.0f;

            // normalize to [0,1]
            float t = glm::clamp(u2 / u2Max, 0.0f, 1.0f);

            // colors
            glm::vec3 cold(0.2f, 0.4f, 1.0f); // blue
            glm::vec3 mid(0.2f, 1.0f, 0.2f); // green
            glm::vec3 hot(1.0f, 0.2f, 0.1f); // red

            glm::vec3 color;

            if (t < 0.5f)
            {
                float localT = t / 0.5f;
                color = glm::mix(cold, mid, localT);
            }
            else
            {
                float localT = (t - 0.5f) / 0.5f;
                color = glm::mix(mid, hot, localT);
            }
                        
            data.instanceColor = color;

            instances.push_back(data);
        }
        // update instance buffer
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, instances.size() * sizeof(InstanceData), instances.data());

        // particle model
        glm::mat4 particleModel = glm::mat4(1.0f);
        particleShader.setMat4("model", particleModel);

        // draw particles
        glBindVertexArray(partVAO);
        glDrawElementsInstanced(GL_TRIANGLES, partIndices.size(), GL_UNSIGNED_INT, 0, instances.size());
        glBindVertexArray(0);


        // ========== cube ===========
        // reset model for cube 
        if (showCube)
        {
            cubeShader.use(); // change to cubeShader 
            glm::mat4 cubeModel = glm::mat4(1.0f);
            cubeModel = glm::translate(cubeModel, glm::vec3(0.0f, 0.0f, 0.0f));  // cube pos
            cubeModel = glm::scale(cubeModel, glm::vec3(xMax, yMax, zMax));
            cubeShader.setMat4("model", cubeModel);
            // send view, projection to cubeShader
            cubeShader.setMat4("projection", projection);
            cubeShader.setMat4("view", view);

            // draw cube 
            glBindVertexArray(cubeVAO);
            glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
            glBindVertexArray(0);
        }        


        // ====== ImGui window =======
        if (showWindow)
        {
            ImGui::Begin("Settings", &showWindow);

            // start - stop simulation
            if (ImGui::Button("Start"))
                startSimulation = true;
            if (ImGui::Button("Stop "))
                startSimulation = false;

            // render cube
            ImGui::Checkbox("Show cube", &showCube);

            // update imGui temp, pressure, energy
            if (windowStep % 20 == 0)
            {
                windowTemp = temperature;
                windowPressure = pressure;
                windowEnergy = Ekin + Epot;
            }
            windowStep++;

            // Energy, temp
            ImGui::Text("Energy : %.1f (lj)", windowEnergy);
            ImGui::Text("Temperature : %.3f (lj)", windowTemp);
            ImGui::Text("Pressure : %.3f (lj)", windowPressure);

            // fps 
            ImGui::Text("FPS : %.1f", ImGui::GetIO().Framerate);

            // close button
            if (ImGui::Button("Close"))
            {
                showWindow = false;
                ImGui_ImplGlfwGL3_Shutdown();
                ImGui::DestroyContext();
                glfwTerminate();
            }

            // close window
            if (showCloseWindow)
            {
                ImGui::Begin("Simulation ended", &showCloseWindow);
                ImGui::Text("The max step has been reached!");
                startSimulation = false;
                
                meltFile.close();   // close file

                ImGui::End();
            }
            ImGui::End();
        }

        ImGui::Render();
        ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }


    glDeleteVertexArrays(1, &partVAO);
    glDeleteBuffers(1, &partVBO);
    glDeleteBuffers(1, &partEBO);
    glDeleteBuffers(1, &instanceVBO);

    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}

// =====================================
// ========= FUNCTIONS =================
// take input from keyboard
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// whenever the window size changed this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

// function for zoom scroll
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    float zoomSensitivity = 3.0f;

    fov -= (float)yoffset * zoomSensitivity;
    if (fov < 20.0f) fov = 20.0f;
    if (fov > 90.0f) fov = 90.0f;
}

// check when right mouse is pressed
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        if (action == GLFW_PRESS)
        {
            rightMousePressed = true;

            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

            glfwGetCursorPos(window, &lastMouseX, &lastMouseY);
        }
        else if (action == GLFW_RELEASE)
        {
            rightMousePressed = false;

            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        }
    }
}

// move camera when right mouse pressed
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    // give priority to imGui when mouse not pressed
    if (ImGui::GetIO().WantCaptureMouse && !rightMousePressed)
        return;
    
    // dont move camera when mouse not pressed
    if (!rightMousePressed)
        return;
    
    float xpos = (float)xposIn;
    float ypos = (float)yposIn;

    float xoffset = xpos - lastMouseX;
    float yoffset = lastMouseY - ypos;

    lastMouseX = xpos;
    lastMouseY = ypos;

    float sensitivity = 0.2f;
    yaw += xoffset * sensitivity;
    pitch += yoffset * sensitivity;

    if (pitch > 89.0f) pitch = 89.0f;
    if (pitch < -89.0f) pitch = -89.0f;
}

// make sphere particle
void generateSphere(float radius, unsigned int sectors, unsigned int stacks, std::vector<particleVertex>& vertices,
                    std::vector<unsigned int>& indices)
{
    const float PI = 3.14159265359f;

    float sectorStep = 2 * PI / sectors;
    float stackStep = PI / stacks;

    for (unsigned int i = 0; i <= stacks; ++i)
    {
        float stackAngle = PI / 2 - i * stackStep;
        float xy = radius * cosf(stackAngle);
        float z = radius * sinf(stackAngle);

        for (unsigned int j = 0; j <= sectors; ++j)
        {
            float sectorAngle = j * sectorStep;

            particleVertex vertex;

            // Position
            vertex.x = xy * cosf(sectorAngle);
            vertex.y = xy * sinf(sectorAngle);
            vertex.z = z;

            // Normal
            vertex.nx = vertex.x / radius;
            vertex.ny = vertex.y / radius;
            vertex.nz = vertex.z / radius;

            vertices.push_back(vertex);
        }
    }

    // Indices
    for (unsigned int i = 0; i < stacks; ++i)
    {
        unsigned int k1 = i * (sectors + 1);
        unsigned int k2 = k1 + sectors + 1;

        for (unsigned int j = 0; j < sectors; ++j)
        {
            if (i != 0)
            {
                indices.push_back(k1 + j);
                indices.push_back(k2 + j);
                indices.push_back(k1 + j + 1);
            }

            if (i != (stacks - 1))
            {
                indices.push_back(k1 + j + 1);
                indices.push_back(k2 + j);
                indices.push_back(k2 + j + 1);
            }
        }
    }
}

// calc random float in [start, finish]
float randomFloat(float start, float finish)
{
    static std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<float> dist(start, finish);
    return dist(gen);
}


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

// particle-particle lj force
glm::vec3 computeLJForce(const glm::vec3& pi, const glm::vec3& pj, float epsilon, float sigma, float& Epot) 
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

// make neighboor list
void buildNeighborList(float rlist)
{
    float rlist2 = rlist * rlist;

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

// calc all lj forces and Epot and Virial 
void computeForces_Epot_Virial(float epsilon, float sigma, float& Epot, float& Virial)
{
    // set forces = 0
    for (auto& p : particles)
        p.force = glm::vec3(0.0f);

    // calc new lj forces
    for (int i = 0; i < N; i++)
    {
        for (int j : neighborList[i])
        {
            if (i >= j) continue;

            glm::vec3 ijForce = computeLJForce(particles[i].position, particles[j].position, epsilon, sigma, Epot);

            particles[i].force += ijForce;
            particles[j].force -= ijForce;


            glm::vec3 rij = particles[i].position - particles[j].position;
            Virial += glm::dot(rij, ijForce);
        }
    }
}

// periodic boundary conditions
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

// check if neighboor list needs rebuild
// IDEA: If any particle has moved (distance > rs / 2), we rebuild the neighborList
bool needRebuild(float skin)
{
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

// calc T (kb = 1, lj units)
float computeTemperature(const float Ekin, const float dof)
{   
    return (2.0f * Ekin) / dof;
}

// thermostat (Berendsen)
void useThermostat(const float dt, const float tau, const float Ttarget, 
                   const float dof, const float Ekin, const float T)
{
    if (T < 1e-8f) return; // avoid division by zero

    // scaling factor
    float lambda = sqrt(1.0f + (dt / tau) * (Ttarget / T - 1.0f));

    // change velocity of particles
    for (auto& p : particles)
        p.velocity *= lambda;
}

// calc pressure (= density * temp + sum(Rij*Fij) / (3V) )
float computePressure(float density, float temp, float virial, float vol)
{
    return density * temp + virial / (3.0f * vol);
}