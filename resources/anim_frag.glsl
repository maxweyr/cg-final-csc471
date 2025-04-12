#version 410 core

in vec3 fragNor;
in vec2 fragTex;
in vec3 fragPos;
in vec3 fragTangent;
in vec3 fragBitangent;

out vec4 color;

// Material structure
struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    vec3 emissive;
    float shininess;
    float opacity;
    bool hasTexture;
};

uniform Material material;
uniform vec3 lightPos;

// Texture samplers
uniform sampler2D texture_diffuse1;
uniform sampler2D texture_specular1;
uniform sampler2D texture_normal1;
uniform sampler2D texture_height1;  // For parallax mapping if implemented
uniform sampler2D texture_emissive1;


void main() {
    // Normalize our vectors
    vec3 normal = normalize(fragNor);
    vec3 lightDir = normalize(lightPos - fragPos);
    vec3 viewDir = normalize(-fragPos);
    
    // Apply normal mapping if available
    if(material.hasTexture && textureSize(texture_normal1, 0).x > 1) {
        // Build TBN matrix
        vec3 T = normalize(fragTangent);
        vec3 B = normalize(fragBitangent);
        vec3 N = normalize(fragNor);
        mat3 TBN = mat3(T, B, N);
        
        // Get normal from normal map and transform to world space
        normal = texture(texture_normal1, fragTex).rgb;
        normal = normal * 2.0 - 1.0;  // Transform from [0,1] to [-1,1]
        normal = normalize(TBN * normal);
    }
    
    // Ambient
    vec3 ambient = material.ambient;
    
    // Get base color from texture or material
    vec3 baseColor;
    if(material.hasTexture && textureSize(texture_diffuse1, 0).x > 1) {
        baseColor = texture(texture_diffuse1, fragTex).rgb;
    } else {
        baseColor = material.diffuse;
    }
    
    // Diffuse
    float diff = max(dot(normal, lightDir), 0.0);
    vec3 diffuse = diff * baseColor;
    
    // Specular (Blinn-Phong)
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);
    
    // Get specular intensity from texture or material
    vec3 specColor;
    if(material.hasTexture && textureSize(texture_specular1, 0).x > 1) {
        specColor = texture(texture_specular1, fragTex).rgb;
    } else {
        specColor = material.specular;
    }
    vec3 specular = spec * specColor;
    
    // Emissive (self-illumination)
    vec3 emissive;
    if(material.hasTexture && textureSize(texture_emissive1, 0).x > 1) {
        emissive = texture(texture_emissive1, fragTex).rgb;
    } else {
        emissive = material.emissive;
    }
    
    // Final color
    vec3 result = ambient + diffuse + specular + emissive;
    
    // Apply opacity
    float alpha = material.opacity;
    if(material.hasTexture && textureSize(texture_diffuse1, 0).x > 1) {
        alpha = texture(texture_diffuse1, fragTex).a * material.opacity;
    }
    
    color = vec4(result, alpha);
}