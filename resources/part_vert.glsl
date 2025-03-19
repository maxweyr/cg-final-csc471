#version 330 core

layout(location = 0) in vec3 vertPos;
layout(location = 1) in vec4 vertColor;

uniform mat4 P;
uniform mat4 M;
uniform mat4 V;

out vec4 partCol;

void main()
{
	mat4 M0 = M;

	M0[0] = vec4(1.0, 0.0, 0.0, 0.0);
	M0[1] = vec4(0.0, 1.0, 0.0, 0.0);
	M0[2] = vec4(0.0, 0.0, 1.0, 0.0);

	gl_Position = P * V * M0 * vec4(vertPos.xyz, 1.0);

    partCol = vertColor;
}
