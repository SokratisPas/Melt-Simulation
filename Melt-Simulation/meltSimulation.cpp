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
#include "headerFiles/particle-system.h"

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
// - Use SCROLL to zoom.
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

// =========== Particle System ========
ParticleSystem partSystem(N);

float particleRadius = 0.3f;    // scale particles (for visuals)

// ==== Instance data struct =====
struct InstanceData {
    glm::vec3 instancePos;
    glm::vec3 instanceColor;
};
std::vector<InstanceData> instances(N);

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
void generateSphere(float radius, unsigned int sectors, unsigned int stacks, std::vector<particleVertex>& vertices,
    std::vector<unsigned int>& indices);

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
    Shader particleShader("..\\Melt-Simulation\\vertexShaders\\meltSimulation.vert",
                          "..\\Melt-Simulation\\fragmentShaders\\meltSimulation.frag");

    Shader cubeShader("..\\Melt-Simulation\\vertexShaders\\meltSimulation2.vert",
                      "..\\Melt-Simulation\\fragmentShaders\\meltSimulation2.frag");

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
    double e = 1.0;             // epsilon (lj units)
    double sigma = 1.0;         // sigma (lj units)

    double r0 = pow(2.0, 1.0 / 6.0) * sigma;
    double a = sqrt(2.0) * r0;  // lattice constant ( = sqrt(2) * 1.122)

    double dt = 0.0005;         // delta time

    float Ttarget = 1.0f;       // target temperature (lj)
    float tau = 0.1f;           // relaxation strength 


    // ==== initialize parameteres ====
    partSystem.initializeParameteres(a);
    
    // remove net momentum
    // it should : sum(vel[i] == 0)
    partSystem.removeNetMomentum();

    // ==== rescale cube =======
    // (rescale cube BEFORE initialize neighborList and computeForces)
    int cells_per_dim = std::ceil(std::cbrt(N / 4.0));
    float L = cells_per_dim * a; // cubes length
    
    partSystem.xMin = partSystem.yMin = partSystem.zMin = 0.0f;
    partSystem.xMax = partSystem.yMax = partSystem.zMax = L;

    // volume, density calculation 
    partSystem.calcVolumeDensity();

    // ==== neighbor list =====
    float rc = 2.5f * sigma;    // r cutoff = 2.5 * 1.0 = 2.5
    float rSkin = 0.3f;         // r skin
    float rlist = rc + rSkin;   // neighbor list radius

    // apply pbc
    for (auto& p : partSystem.particles)
    {
        partSystem.applyPBC(p.position);
    }

    // initalize neighbour list and forces
    partSystem.buildNeighborList(rlist);
    partSystem.computeForces_Epot_Virial(e, sigma);

    // initiallize acceleration
    partSystem.initiallizeAcc();

    // change target for orbital camera
    target = glm::vec3(partSystem.xMax / 2);
    
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
            partSystem.resetEkinEpotVirial();

            // update pos wth velocity Verlet method
            partSystem.updatePos(dt);

            // check if rebuild neighboorList
            if (partSystem.needRebuild(rSkin))
            {
                partSystem.buildNeighborList(rlist);
            }

            // update forces and E potential and Virial
            partSystem.computeForces_Epot_Virial(e, sigma);

            // acceleration update
            partSystem.updateAcc();

            // vel and Ekin update
            partSystem.updateVelEkin(dt);

            // Temperature
            partSystem.computeTemperature();

            // Pressure
            partSystem.computePressure();

            // thermostat
            if (USE_THERMOSTAT)
            {
                partSystem.useThermostat(dt, tau, Ttarget);
            }

            // ======= WRITE FILE ========            
            if (step % 10 == 0)
            {
                float Etot = partSystem.Ekin + partSystem.Epot;

                meltFile << step << " "
                    << simTime << " "
                    << partSystem.Ekin << " "
                    << partSystem.Epot << " "
                    << Etot << " "
                    << partSystem.temperature << " "
                    << partSystem.pressure << " "
                    << deltaTime << "\n";
            }
            step++;
            simTime += dt;

            if (step == stepMax + 1)
                showCloseWindow = true;
        }


        // update instance pos and color
        instances.clear();
        for (auto& part : partSystem.particles)
        {
            InstanceData data;
            // position
            data.instancePos = part.position;

            // velocity^2
            float u2 = glm::dot(part.velocity, part.velocity);

            // average v^2 from temperature
            float u2Mean = 3.0f * partSystem.temperature;

            // choose upper visualization limit
            float u2Max = 4.0f * u2Mean;

            // avoid divide by zero
            if (u2Max < 1e-6f)
                u2Max = 1.0f;

            // normalize to [0,1]
            float t = glm::clamp(u2 / u2Max, 0.0f, 1.0f);

            // colors
            glm::vec3 cold(0.2f, 0.4f, 1.0f);   // blue
            glm::vec3 mid(0.2f, 1.0f, 0.2f);    // green
            glm::vec3 hot(1.0f, 0.2f, 0.1f);    // red

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
            cubeModel = glm::scale(cubeModel, glm::vec3(partSystem.xMax, partSystem.yMax, partSystem.zMax));
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
                windowTemp = partSystem.temperature;
                windowPressure = partSystem.pressure;
                windowEnergy = partSystem.Ekin + partSystem.Epot;
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