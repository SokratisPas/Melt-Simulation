// meltSimulation.cpp
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
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
#include "headerFiles/orbit-Camera.h"
#include "headerFiles/POSCAR-parser.h"


// ==== particle vertex struct =====
struct particleVertex {
    float x, y, z;      // pos
    float nx, ny, nz;   // normals
};

// =========== Particle System ========
ParticleSystem partSystem(1);

float particleRadius = 0.3f;    // scale particles (for visuals)


// Starting screen settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 750;

// ==== orbit camera =====
OrbitCamera cam;

bool rightMousePressed = false;

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
    // POSCAR FILE
    POSCAR poscar;
    if (!poscar.load("poscar-input1.txt"))
        return 1;

    int N = poscar.atoms.size();
    std::cout << N;

    for (size_t i = 0; i < poscar.atoms.size(); ++i)
    {
        std::cout
            << i << " "
            << poscar.atoms[i].element << " "
            << poscar.atoms[i].x << " "
            << poscar.atoms[i].y << " "
            << poscar.atoms[i].z << '\n';
    }

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
    bool showAxis = false;

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

    Shader axisShader("..\\Melt-Simulation\\vertexShaders\\meltSimulation2.vert",
                      "..\\Melt-Simulation\\fragmentShaders\\meltSimulation2.frag");

    // =================================
    // =========== MODEL DATA ==========
    //==================================
    // =========== Particle =============
    std::vector<particleVertex> partVerts;
    std::vector<unsigned int> partIndices;
    unsigned int partSectors = 22;
    unsigned int partStacks = 16;
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


    // =========== Axis ============
    float axisVertices[] =
    {
        // X axis
        0.0f, 0.0f, 0.0f,
        5.0f, 0.0f, 0.0f,

        // Y axis
        0.0f, 0.0f, 0.0f,
        0.0f, 5.0f, 0.0f,

        // Z axis
        0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 5.0f
    };

    unsigned int axisVAO, axisVBO;

    glGenVertexArrays(1, &axisVAO);
    glGenBuffers(1, &axisVBO);

    glBindVertexArray(axisVAO);

    glBindBuffer(GL_ARRAY_BUFFER, axisVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(axisVertices), axisVertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);


    // ======================================
    // =========== PHYSICS DATA =============

    // ==== initialize parameteres ====
    partSystem.particles.resize(N);


    for (int i = 0; i < N; i++)
    {
        // initial pos
        partSystem.particles[i].position = glm::vec3(poscar.atoms[i].x, poscar.atoms[i].y, poscar.atoms[i].z);
    }
    

    // change target for orbital camera
    glm::vec3 center(0.0f);

    for (const auto& atom : poscar.atoms)
    {
        center += glm::vec3(atom.x, atom.y, atom.z);
    }

    center /= static_cast<float>(poscar.atoms.size());

    cam.target = center;
    

    // enable depth test
    glEnable(GL_DEPTH_TEST); 

    // face culling
    glEnable(GL_CULL_FACE);

    // anti aliasing
    glEnable(GL_MULTISAMPLE);

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
        glm::mat4 projection = cam.calcProjection(fullWidth, fullHeight);
        
        // view for particleShader    
        glm::mat4 view = cam.calcView();

        particleShader.setMat4("projection", projection);
        particleShader.setMat4("view", view);

        for (const auto& atom : poscar.atoms)
        {
            glm::mat4 model(1.0f);

            // particle model
            glm::mat4 particleModel = glm::mat4(1.0f);
            particleModel = glm::translate(particleModel, glm::vec3(atom.x, atom.y, atom.z));
            particleModel = glm::scale(particleModel, glm::vec3(particleRadius));
            particleShader.setMat4("model", particleModel);

            // color
            glm::vec3 color(1.0f);
            if (atom.element == "Ti")
                color = glm::vec3(1.0f, 0.0f, 0.0f);
            else if (atom.element == "Sr")
                color = glm::vec3(0.0f, 1.0f, 0.0f);

            particleShader.setVec3("color", color);

            // draw particles
            glBindVertexArray(partVAO);
            glDrawElements(GL_TRIANGLES, partIndices.size(), GL_UNSIGNED_INT, 0);
            glBindVertexArray(0);
        }

        // ========== Axis ===========
        if (showAxis)
        {
            axisShader.use();
            axisShader.setMat4("projection", projection);
            axisShader.setMat4("view", view);

            glm::mat4 axisModel = glm::mat4(1.0f);
            axisModel = glm::translate(axisModel, center);
            axisModel = glm::scale(axisModel, glm::vec3(10.0f));
            axisShader.setMat4("model", axisModel);

            glBindVertexArray(axisVAO);
            glDrawArrays(GL_LINES, 0, 6);
            glBindVertexArray(0);
        }        


        // ====== ImGui window =======
        if (showWindow)
        {
            ImGui::Begin("Settings", &showWindow);

            // show Axis
            ImGui::Checkbox("Show Axis", &showAxis);

            // change partile scale 
            ImGui::SliderFloat("Radius", &particleRadius, 0.1f, 3.0f);

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

    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}

// =====================================
// ========= FUNCTIONS =================
// 
// -----------------------------------------
// take input from keyboard
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// -----------------------------------------
// whenever the window size changed this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

// -----------------------------------------
// function for zoom scroll
// here the xoffset is not needed, but the "glfwSetScrollCallback" function needs it 
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    cam.zoomScroll(yoffset);
}

// -----------------------------------------
// check when right mouse is pressed
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        if (action == GLFW_PRESS)           // right mouse pressed
        {
            rightMousePressed = true;

            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED); // hide cursor

            glfwGetCursorPos(window, &cam.lastMouseX, &cam.lastMouseY); // update last mouse pos
        }
        else if (action == GLFW_RELEASE)    // right mouse released
        {
            rightMousePressed = false;

            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL); // show again cursor
        }
    }
}

// -----------------------------------------
// move camera when right mouse pressed
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    // give priority to imGui when mouse not pressed
    if (ImGui::GetIO().WantCaptureMouse && !rightMousePressed)
        return;
    
    // dont move camera when mouse not pressed
    if (!rightMousePressed)
        return;
    
    // move camera
    cam.moveCamera(xposIn, yposIn);
}

// -----------------------------------------
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