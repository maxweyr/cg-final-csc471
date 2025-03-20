#version 330 core

in vec3 fragNor;
in vec2 fragTex;

out vec4 color;

uniform vec3 MatAmb;
uniform vec3 MatDif;
uniform vec3 MatSpec;
uniform float MatShine;
uniform vec3 lightPos;
uniform int hasTexture;
uniform sampler2D textureSampler;

void main() {
    vec3 normal = normalize(fragNor);
    vec3 lightDir = normalize(lightPos);
    
    // Ambient
    vec3 ambient = MatAmb;
    
    // Get base color from texture or material
    vec3 baseColor;
    if(hasTexture == 1) {
        baseColor = texture(textureSampler, fragTex).rgb;
    } else {
        baseColor = MatDif;
    }
    
    // Diffuse
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = diff * baseColor;
    
    // Specular (Blinn-Phong)
    vec3 viewDir = vec3(0.0, 0.0, 1.0); // Simplified view direction
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), MatShine);
    vec3 specular = spec * MatSpec;
    
    // Final color
    vec3 result = ambient + diffuse + specular;
    color = vec4(result, 1.0);
}