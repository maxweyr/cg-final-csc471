#version 410 core
layout(location = 0) in vec3 vertPos;
layout(location = 1) in vec3 vertNor;
uniform mat4 P;
uniform mat4 V;
uniform mat4 M;
out vec3 fragNor;
void main() {
    // Simple transformation
    vec4 pos = vec4(vertPos, 1.0);
    
    // Transform normal
    mat3 normalMatrix = transpose(inverse(mat3(M)));
    fragNor = normalMatrix * vertNor;
    
    // Final position
    gl_Position = P * V * M * pos;
}