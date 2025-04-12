#version 410 core

layout(location = 0) in vec3 vertPos;
layout(location = 1) in vec3 vertNor;
layout(location = 2) in vec2 vertTex;
layout(location = 3) in vec3 vertTan;    // New: tangent attribute
layout(location = 4) in vec3 vertBitan;  // New: bitangent attribute
layout(location = 5) in ivec4 boneIDs;   // Changed layout location
layout(location = 6) in vec4 weights;    // Changed layout location

uniform mat4 P;
uniform mat4 V;
uniform mat4 M;

#define MAX_BONES 100
uniform mat4 boneTransforms[MAX_BONES];

out vec3 fragNor;
out vec2 fragTex;
out vec3 fragPos;       // New: fragment position for lighting calculations
out vec3 fragTangent;   // New: tangent vector for normal mapping
out vec3 fragBitangent; // New: bitangent vector for normal mapping

void main() {
    // Apply bone transformations
    mat4 boneTransform = mat4(0.0);
    
    for(int i = 0; i < 4; i++) {
        if(boneIDs[i] != -1) {
            boneTransform += weights[i] * boneTransforms[boneIDs[i]];
        }
    }
    
    // If no bones affect this vertex or failed to load, use identity
    if(boneTransform == mat4(0.0)) {
        boneTransform = mat4(1.0);
    }
    
    // Transform position with bone transform
    vec4 pos = boneTransform * vec4(vertPos, 1.0);
    
    // Transform normal with bone transform (normal matrix)
    mat3 normalMatrix = transpose(inverse(mat3(boneTransform * M)));
    fragNor = normalMatrix * vertNor;
    
    // Transform tangent and bitangent
    fragTangent = normalMatrix * vertTan;
    fragBitangent = normalMatrix * vertBitan;
    
    // Pass through texture coordinates
    fragTex = vertTex;
    
    // Calculate fragment position in world space for lighting
    fragPos = vec3(M * pos);
    
    // Final position
    gl_Position = P * V * M * pos;
}