#version 330 core

// Texture sampler
uniform sampler2D Texture;
uniform float MatShine;
uniform int flip;   // 0 = regular normals, 1 = flipped normals
uniform int ground; // 0 = no hole, # = hole of radius #
uniform int basin;  // 0 = normal drawing, 1 = only draw below ground plane

in vec2 vTexCoord;
in vec3 fragNor;
in vec3 lightDir;
in vec3 EPos;
in vec3 fragWorldPos;

out vec4 Outcolor;

void main() {
    // hole test
    if (ground > 0) {
        // Calculate distance from center in world space
        float dist = length(vec2(fragWorldPos.x, fragWorldPos.z));
    
        // Discard fragment if inside the hole
        if (dist < ground) {
            discard;
        }
    }

    // Basin test - only draw below ground plane
    if (basin > 0 && fragWorldPos.y > -1.25) {
        discard;
    }

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