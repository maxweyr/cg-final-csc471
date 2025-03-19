#version 330 core

layout(location = 0) in vec4 vertPos;
layout(location = 1) in vec3 vertNor;
layout(location = 2) in vec2 vertTex;

uniform mat4 P;
uniform mat4 V;
uniform mat4 M;

uniform vec3 lightPos;
uniform int useFixedLight;
uniform vec3 fixedLightPos;

uniform float uTime;         // time in seconds
uniform float windStrength;  // how far vertices are displaced
uniform float windFrequency; // frequency of the sine wave
uniform bool applyWind;      // flag to apply wind effect

// Outputs
out vec3 fragNor;
out vec3 lightDir;
out vec3 EPos;

void main()
{
    vec4 pos = vertPos;
    if (applyWind) {
        float displacement = sin(uTime + vertPos.y * windFrequency) * windStrength;
        pos.x += displacement;
        pos.z += displacement;
    }
    gl_Position = P * V * M * pos;
    
    fragNor = (M * vec4(vertNor, 0.0)).xyz;

    if (useFixedLight == 1) {
		lightDir = fixedLightPos - (M * vertPos).xyz;
	} else {
        lightDir = lightPos - (M * vertPos).xyz;
    }
    

    EPos = (M * vertPos).xyz;
}