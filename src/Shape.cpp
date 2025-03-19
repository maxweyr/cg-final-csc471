
#include "Shape.h"
#include <iostream>
#include <cassert>
#include <limits>
#include "GLSL.h"
#include "Program.h"

using namespace std;


// copy the data from the shape to this object
void Shape::createShape(tinyobj::shape_t & shape)
{
	posBuf = shape.mesh.positions;
	norBuf = shape.mesh.normals;
	texBuf = shape.mesh.texcoords;
	eleBuf = shape.mesh.indices;

	if (texBuf.empty()) {
		cout << "Warning: No texture coordinates found for shape!" << endl;
	}
}

void Shape::measure()
{
	float minX, minY, minZ;
	float maxX, maxY, maxZ;

	minX = minY = minZ = (std::numeric_limits<float>::max)();
	maxX = maxY = maxZ = -(std::numeric_limits<float>::max)();

	//Go through all vertices to determine min and max of each dimension
	for (size_t v = 0; v < posBuf.size() / 3; v++)
	{
		if (posBuf[3*v+0] < minX) minX = posBuf[3 * v + 0];
		if (posBuf[3*v+0] > maxX) maxX = posBuf[3 * v + 0];

		if (posBuf[3*v+1] < minY) minY = posBuf[3 * v + 1];
		if (posBuf[3*v+1] > maxY) maxY = posBuf[3 * v + 1];

		if (posBuf[3*v+2] < minZ) minZ = posBuf[3 * v + 2];
		if (posBuf[3*v+2] > maxZ) maxZ = posBuf[3 * v + 2];
	}

	min.x = minX;
	min.y = minY;
	min.z = minZ;
	max.x = maxX;
	max.y = maxY;
	max.z = maxZ;
}

void Shape::init()
{
	// Initialize the vertex array object
	CHECKED_GL_CALL(glGenVertexArrays(1, &vaoID));
	CHECKED_GL_CALL(glBindVertexArray(vaoID));

	// Send the position array to the GPU
	CHECKED_GL_CALL(glGenBuffers(1, &posBufID));
	CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, posBufID));
	CHECKED_GL_CALL(glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STATIC_DRAW));

	// Send the normal array to the GPU
	if (norBuf.empty()) {	// check for normals
		cerr << "Warning: Found no normals! Generating..." << endl;

		norBuf.resize(posBuf.size(), 0.0f);	// Resize the normal buffer to match the position buffer.

		// compute face normals (assuming triangles)
		for (size_t i = 0; i < eleBuf.size(); i += 3) {
			// Retrieve indices for the triangle.
			int idx0 = eleBuf[i];
			int idx1 = eleBuf[i + 1];
			int idx2 = eleBuf[i + 2];

			// Get the three vertex positions.
			glm::vec3 p0(posBuf[3 * idx0], posBuf[3 * idx0 + 1], posBuf[3 * idx0 + 2]);
			glm::vec3 p1(posBuf[3 * idx1], posBuf[3 * idx1 + 1], posBuf[3 * idx1 + 2]);
			glm::vec3 p2(posBuf[3 * idx2], posBuf[3 * idx2 + 1], posBuf[3 * idx2 + 2]);

			// compute face normal
			glm::vec3 normal = glm::normalize(glm::cross(p1 - p0, p2 - p0));

			// accumulate normals into each vertex's normal
			norBuf[3 * idx0] += normal.x;
			norBuf[3 * idx0 + 1] += normal.y;
			norBuf[3 * idx0 + 2] += normal.z;

			norBuf[3 * idx1] += normal.x;
			norBuf[3 * idx1 + 1] += normal.y;
			norBuf[3 * idx1 + 2] += normal.z;

			norBuf[3 * idx2] += normal.x;
			norBuf[3 * idx2 + 1] += normal.y;
			norBuf[3 * idx2 + 2] += normal.z;
		}

		// normalize all accumulated normals for each vertex
		for (size_t i = 0; i < posBuf.size() / 3; i++)
		{
			glm::vec3 n(norBuf[3 * i], norBuf[3 * i + 1], norBuf[3 * i + 2]);
			n = glm::normalize(n);
			norBuf[3 * i] = n.x;
			norBuf[3 * i + 1] = n.y;
			norBuf[3 * i + 2] = n.z;
		}

		norBufID = 1;
	}
	CHECKED_GL_CALL(glGenBuffers(1, &norBufID));
	CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, norBufID));
	CHECKED_GL_CALL(glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_STATIC_DRAW));

	// Send the texture array to the GPU
	if (texBuf.empty())
	{
		texBufID = 0;
		cout << "Warning: no textures!" << endl;
	}
	else
	{
		CHECKED_GL_CALL(glGenBuffers(1, &texBufID));
		CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, texBufID));
		CHECKED_GL_CALL(glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW));
	}

	// Send the element array to the GPU
	CHECKED_GL_CALL(glGenBuffers(1, &eleBufID));
	CHECKED_GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID));
	CHECKED_GL_CALL(glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size()*sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW));

	// Unbind the arrays
	CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
	CHECKED_GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
}

void Shape::draw(const shared_ptr<Program> prog) const
{
	int h_pos, h_nor, h_tex;
	h_pos = h_nor = h_tex = -1;

	CHECKED_GL_CALL(glBindVertexArray(vaoID));

	// Bind position buffer
	h_pos = prog->getAttribute("vertPos");
	GLSL::enableVertexAttribArray(h_pos);
	CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, posBufID));
	CHECKED_GL_CALL(glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0));

	// Bind normal buffer
	h_nor = prog->getAttribute("vertNor");
	if (h_nor != -1 && norBufID != 0)
	{
		GLSL::enableVertexAttribArray(h_nor);
		CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, norBufID));
		CHECKED_GL_CALL(glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0));
	}

	// check if texture coordinates are available AND if the shader has the vertTex attribute
	h_tex = prog->getAttribute("vertTex");
	if (h_tex != -1 && texBufID != 0)
	{
		GLSL::enableVertexAttribArray(h_tex);
		CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, texBufID));
		CHECKED_GL_CALL(glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void *)0));
	}

	// Bind element buffer
	CHECKED_GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID));

	// Draw
	CHECKED_GL_CALL(glDrawElements(GL_TRIANGLES, (int)eleBuf.size(), GL_UNSIGNED_INT, (const void *)0));

	// Disable and unbind
	if (h_tex != -1 && texBufID != 0)
	{
		GLSL::disableVertexAttribArray(h_tex);
	}
	if (h_nor != -1 && norBufID != 0)
	{
		GLSL::disableVertexAttribArray(h_nor);
	}
	GLSL::disableVertexAttribArray(h_pos);
	CHECKED_GL_CALL(glBindBuffer(GL_ARRAY_BUFFER, 0));
	CHECKED_GL_CALL(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
}
