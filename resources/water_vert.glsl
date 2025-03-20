#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 P;
uniform mat4 V;
uniform mat4 M;
uniform float uTime;

out vec3 FragPos;

void main() {
    // Calculate displacement based on time and position
    vec3 pos = aPos;
    
    // Apply vertical displacement for actual geometric waves
    float wave1 = sin(pos.x * 2.0 + uTime * 1.2) * 0.1;
    float wave2 = sin(pos.z * 3.0 + uTime * 0.8) * 0.1;
    float wave3 = sin((pos.x + pos.z) * 1.5 + uTime * 1.5) * 0.05;
    
    // Apply the displacement
    pos.y += wave1 + wave2 + wave3;
    
    gl_Position = P * V * M * vec4(pos, 1.0);
    FragPos = vec3(M * vec4(pos, 1.0));
}