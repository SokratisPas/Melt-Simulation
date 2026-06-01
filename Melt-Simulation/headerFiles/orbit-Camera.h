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

	float distanceToTarget = 50.0f;
	glm::vec3 target = glm::vec3(0.5f, 0.5f, 0.5f); // target where camera looks

	float cameraSensitivity = 0.2f;
	float zoomSensitivity = 3.0f;

	// (move camera some offset in yaw and pitch planes, so we need the prev positions)
	double lastMouseX = 0.0f;	// last mouse x pos
	double lastMouseY = 0.0f;	// last mouse y pos

	bool showOrtho = false;		// show orthocanonical projection

	float orthoSize = 30.0f;	// size for orth projection

	// ================================
	// ======== functions =============
	// 
	// ---------------------------------------
	// projection
	glm::mat4 calcProjection(int width, int height)
	{
		float aspectRatio = (float)width / (float)height;

		if (showOrtho)
			return glm::ortho(-orthoSize * aspectRatio, orthoSize * aspectRatio, -orthoSize, orthoSize, 0.1f, 200.0f);		

		return glm::perspective(glm::radians(fov), aspectRatio, 0.1f, 200.0f);
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

	// --------------------------------------
	// zoom with keyboard
	void zoomWith_Z()
	{
		if (showOrtho)	// ortho zoom
		{
			orthoSize -= 2.0f;
			orthoSize = glm::clamp(orthoSize, 5.0f, 100.0f);
		}
		else {			// perspective zoom
			distanceToTarget -= 0.5f;
			distanceToTarget = glm::clamp(distanceToTarget, 1.0f, 100.0f);
		}
	}

	// --------------------------------------
	// unzoom with keyboard
	void unzoomWith_X()
	{
		if (showOrtho)	// ortho unzoom
		{
			orthoSize += 2.0f;
			orthoSize = glm::clamp(orthoSize, 5.0f, 100.0f);
		}
		else {			// perspective unzoom
			distanceToTarget += 0.5f;
			distanceToTarget = glm::clamp(distanceToTarget, 1.0f, 100.0f);
		}
	}
};