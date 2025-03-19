#version 330 core

in vec3 fragNor;		//interpolated normal vector in camera space
in vec3 lightDir;		//interpolated light vector in camera space
in vec3 EPos;			//position of the vertex in camera space

out vec4 color;

uniform vec3 lightPos;  // Light position in world space
uniform vec3 MatAmb;    // Ambient color
uniform vec3 MatDif;    // Diffuse color
uniform vec3 MatSpec;   // Specular color
uniform float MatShine; // Shininess
uniform vec3 lightIntensity; // Light intensity
uniform int flip; // Uniform to control normal flip
uniform vec3 ambientLight;     // Global ambient light intensity
uniform vec4 starColor;      // Star color and alpha
uniform bool ignoreLight;

void main() {
	vec3 normal = normalize(fragNor);

	// Apply normal flipping ONLY if flip == 1
    if (flip == 1) {
        normal = -normal;
    }

	if (ignoreLight) {
        // for shadows, just use a constant black color with alpha
        color = vec4(0.0, 0.0, 0.0, 0.5); // alpha effects shadow intensity
    } else {
		vec3 light = normalize(lightDir);

		// Ambient component for global illumination
		vec3 ambient = ambientLight * MatAmb;

		//compute view direction
		vec3 view = normalize(-EPos);

		//compute diffuse component
		float dc = max(dot(normal, light), 0.0);
		vec3 diffuse = dc * MatDif * lightIntensity;

		//compute halfway vector
		vec3 h = normalize(light + view);

		//compute specular component
		float sc = pow(max(dot(normal, h), 0.0), MatShine);
		vec3 specular = sc * MatSpec * lightIntensity;

		//compute final color
		//compute with ambient here so that the global light does not effect shadows
		vec3 finalColor = ambient + diffuse + specular;

		//------------------------ set color -------------------------------
		// check if starColor is being used (alpha < 1.0)
		if (starColor.a < 1.0) {
			color = vec4(finalColor, starColor.a);
		} else {
			color = vec4(finalColor, 1.0);
		}
	}
}
