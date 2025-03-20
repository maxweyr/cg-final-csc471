#version 330 core

layout(location = 0) in vec3 vertPos;
layout(location = 1) in vec3 vertNor;
layout(location = 2) in vec2 vertTex;
layout(location = 3) in ivec4 boneIDs;
layout(location = 4) in vec4 weights;

uniform mat4 P;
uniform mat4 V;
uniform mat4 M;
#define MAX_BONES 100
uniform mat4 boneTransforms[MAX_BONES];

out vec3 fragNor;
out vec2 fragTex;

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
    mat3 normalMatrix = transpose(inverse(mat3(boneTransform)));
    fragNor = normalMatrix * vertNor;
    
    // Pass through texture coordinates
    fragTex = vertTex;
    
    // Final position
    gl_Position = P * V * M * pos;
}