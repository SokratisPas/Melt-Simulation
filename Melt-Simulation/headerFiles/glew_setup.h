#ifndef CAMERA_H
#define CAMERA_H

#include<GL/glew.h>

#include<cstdio>
#include<cstdlib>

void glew_setup()
{
    glewExperimental = true;
    if (glewInit() != GLEW_OK)
    {
        printf("Failed to initialize glew. Exiting...");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
}

#endif