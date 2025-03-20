#version 330 core

out vec4 FragColor;
in vec3 FragPos;

uniform float uTime;
uniform vec3 lightPos;
uniform vec3 viewPos;

void main() {
    // Add ripple effect based on position and time
    float ripple1 = sin(FragPos.x * 3.0 + uTime * 1.5) * 0.5 + 0.5;
    float ripple2 = sin(FragPos.z * 2.0 + uTime * 1.0) * 0.5 + 0.5;
    float ripple3 = sin((FragPos.x + FragPos.z) * 2.5 + uTime * 2.0) * 0.5 + 0.5;
    
    float rippleEffect = (ripple1 + ripple2 + ripple3) / 3.0;
    
    // Lighten color based on ripple
    vec3 deepColor = vec3(0.0, 0.2, 0.4);  // Darker blue for troughs
    vec3 shallowColor = vec3(0.1, 0.5, 0.7); // Lighter blue for peaks
    vec3 finalColor = mix(deepColor, shallowColor, rippleEffect);
    
    // Calculate light direction
    vec3 lightDir = normalize(lightPos - FragPos);
    
    // Simple lighting calculation
    float diff = max(dot(vec3(0, 1, 0), lightDir), 0.0);
    vec3 diffuse = diff * vec3(1.0, 1.0, 0.9); // Moonlight color
    
    // View direction
    vec3 viewDir = normalize(viewPos - FragPos);
    
    // Calculate a normal based on ripple gradients
    vec3 rippleNormal = normalize(vec3(
        -(ripple1 * 2.0 - 1.0) * 0.5,
        1.0,
        -(ripple2 * 2.0 - 1.0) * 0.5
    ));
    
    // Add diagonal ripple influence to both X and Z
    vec3 diagonalNormal = normalize(vec3(
        -(ripple3 * 2.0 - 1.0) * 0.3,
        1.0,
        -(ripple3 * 2.0 - 1.0) * 0.3
    ));
    
    // Blend the two normals
    rippleNormal = normalize(rippleNormal + diagonalNormal);
    
    // Use the ripple normal for better specular highlights
    vec3 reflectDir = reflect(-lightDir, rippleNormal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32); // Reduced from 64 to 32 for wider highlights
    
    // Add secondary light contribution from the opposite direction for more consistent highlights
    vec3 secondaryLightDir = normalize(vec3(-lightPos.x, lightPos.y, -lightPos.z) - FragPos);
    vec3 secondaryReflectDir = reflect(-secondaryLightDir, rippleNormal);
    float secondarySpec = pow(max(dot(viewDir, secondaryReflectDir), 0.0), 32);
    
    // Combine the specular contributions
    float combinedSpec = max(spec, secondarySpec * 0.3); // Secondary light is dimmer
    vec3 specular = vec3(1.0) * combinedSpec * 1.0; // Increased from 0.8 to 1.0 for brighter specular
    
    // Enhanced Fresnel effect for more visible reflections
    float fresnelTerm = 0.5 + 0.5 * pow(1.0 - max(dot(viewDir, rippleNormal), 0.0), 3.0);
    
    // Add moonlight reflection based on angle
    vec3 moonReflection = vec3(1.0, 0.98, 0.9) * fresnelTerm * combinedSpec * 0.8;
    finalColor += moonReflection;
    
    // Add a subtle ambient reflection that's always visible regardless of moon position
    vec3 ambientReflection = vec3(0.1, 0.1, 0.15) * fresnelTerm;
    finalColor += ambientReflection;
    
    // Final color with lighting
    finalColor = finalColor * (diffuse + vec3(0.3)) + specular;
    
    // Add some transparency
    FragColor = vec4(finalColor, 0.8);
}