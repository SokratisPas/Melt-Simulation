#pragma once

#include <GL/glew.h>

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtc/type_ptr.hpp>

class OrbitCamera {
public: 
	float yaw = -90.0f;
	float pitch = 0.0f;

	float fov = 45.0f; 

	float distanceToTarget = 20.0f;
	glm::vec3 target = glm::vec3(0.5f, 0.5f, 0.5f); // target where camera looks

	float cameraSensitivity = 0.2f;
	float zoomSensitivity = 3.0f;

	// (move camera some offset in yaw and pitch planes, so we need the prev positions)
	double lastMouseX = 0.0f;	// last mouse x pos
	double lastMouseY = 0.0f;	// last mouse y pos

	// ================================
	// ======== functions =============
	// 
	// ---------------------------------------
	// projection
	glm::mat4 calcProjection(int width, int height)
	{
		return glm::perspective(glm::radians(fov), (float)width / (float)height, 0.1f, 100.0f);
	}

	// ---------------------------------------
	// view
	glm::mat4 calcView()
	{
		glm::vec3 direction;
		direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
		direction.y = sin(glm::radians(pitch));
		direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));

		glm::vec3 cameraPos = target - direction * distanceToTarget;

		return glm::lookAt(cameraPos, target, glm::vec3(0.0f, 1.0f, 0.0f));
	}

	// ---------------------------------------
	// scroll function
	void zoomScroll(double yoffset)
	{
		fov -= (float)yoffset * zoomSensitivity;
		if (fov < 20.0f) fov = 20.0f;
		if (fov > 90.0f) fov = 90.0f;
	}

	// ---------------------------------------
	// move camera
	void moveCamera(double xposIn, double yposIn)
	{
		float xpos = (float)xposIn;
		float ypos = (float)yposIn;

		float xoffset = xpos - lastMouseX;
		float yoffset = lastMouseY - ypos;

		lastMouseX = xpos;
		lastMouseY = ypos;
		
		yaw += xoffset * cameraSensitivity;
		pitch += yoffset * cameraSensitivity;

		if (pitch > 89.0f) pitch = 89.0f;
		if (pitch < -89.0f) pitch = -89.0f;
	}
};