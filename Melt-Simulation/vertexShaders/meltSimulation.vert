#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;


uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

uniform vec3 color;

out vec3 partColor;
out vec3 Normal;
out vec3 FragPos;

void main()
{
    // color
    partColor = color;

    // normals
    FragPos = vec3(model * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(model))) * aNormal;

    // position 
    gl_Position = projection * view * model * vec4(aPos, 1.0);   
}