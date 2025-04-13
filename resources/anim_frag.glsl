#version 410 core
in vec3 fragNor;
out vec4 color;
uniform vec3 solidColor;
void main() {
    // Just use a solid color for debug
    color = vec4(solidColor, 1.0);
}