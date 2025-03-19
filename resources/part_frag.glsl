#version 330 core

uniform sampler2D alphaTexture;

in vec4 partCol;

out vec4 outColor;

void main()
{
	// Sample the alpha texture using built-in point coordinates
	float alpha = texture(alphaTexture, gl_PointCoord).r;

	// Multiply the vertex alpha by the texture's alpha.
    outColor = vec4(partCol.rgb, partCol.a * alpha);
	//outColor = vec4(partCol, alpha);
}
