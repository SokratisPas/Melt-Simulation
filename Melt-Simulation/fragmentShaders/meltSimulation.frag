#version 330 core

in vec3 partColor;
in vec3 Normal;
in vec3 FragPos;

out vec4 FragColor;

void main()
{
    vec3 lightDir = normalize(vec3(0.0, 1.0, 0.0)); // light is from top

    float diff = max(dot(normalize(Normal), lightDir), 0.0);

    float ambient = 0.2f; // ambient light

    vec3 color = partColor * (diff + ambient);

    FragColor = vec4(color, 1.0);
}