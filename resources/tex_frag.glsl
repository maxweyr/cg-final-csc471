#version 330 core

// Texture sampler
uniform sampler2D Texture;
uniform float MatShine;
uniform int flip; // Uniform to control normal flip

in vec2 vTexCoord;  // Correctly passing texture coordinates
in vec3 fragNor;
in vec3 lightDir;
in vec3 EPos;

out vec4 Outcolor;

void main() {
    vec4 texColor = texture(Texture, vTexCoord); // Sample texture

    vec3 normal = normalize(fragNor);
    // Apply normal flipping ONLY if flip == 1
    if (flip == 1) {
        normal = -normal;
    }

    vec3 light = normalize(lightDir);
    // Compute diffuse coefficient
    float dC = max(0.0, dot(normal, light));

    // Compute specular coefficient (Blinn-Phong)
    vec3 viewDir = normalize(-EPos);
    vec3 halfV = normalize(viewDir + light);
    float sC = pow(max(dot(normalize(halfV), normal), 0), MatShine);

    // Apply shading to texture color
    vec3 finalColor = texColor.rgb * dC + texColor.rgb * sC;

    // Output final shaded color
    Outcolor = vec4(finalColor, texColor.a);
}