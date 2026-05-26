#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec3 InstancePos;
layout (location = 3) in vec3 InstanceColor;

out vec3 partColor;
out vec3 Normal;
out vec3 FragPos;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

void main()
{
    // color
    partColor = InstanceColor;

    // normals
    FragPos = vec3(model * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(model))) * aNormal;

    // position    
    vec3 pos = aPos + InstancePos;
    gl_Position = projection * view * model * vec4(pos, 1.0);   
}