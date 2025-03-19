#define GLM_ENABLE_EXPERIMENTAL
#include <iostream>
#include <algorithm>
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/quaternion.hpp>
#include "particleSys.h"
#include "GLSL.h"

using namespace std;

particleSys::particleSys(vec3 source) {
	numP = 300;	
	t = 0.0f;
	h = 0.01f;
	g = vec3(0.0f, -0.098, 0.0f);
	start = vec3(0);
	theCamera = glm::mat4(1.0);
}

void particleSys::gpuSetup() {

    std::cout << "Initializing particle system with " << numP << " particles" << std::endl;
    std::cout << "Emitter position: " << emitterPos.x << ", " << emitterPos.y << ", " << emitterPos.z << std::endl;

    // Initialize particles and arrays
    for (int i = 0; i < numP; i++) {
        points[i * 3 + 0] = emitterPos.x;
        points[i * 3 + 1] = emitterPos.y;
        points[i * 3 + 2] = emitterPos.z;

        auto particle = std::make_shared<Particle>(emitterPos);
        particles.push_back(particle);
        particle->load(emitterPos);

        // Set initial colors
        vec4 col = particle->getColor();
        pointColors[i * 4 + 0] = col.r;
        pointColors[i * 4 + 1] = col.g;
        pointColors[i * 4 + 2] = col.b;
        pointColors[i * 4 + 3] = col.a;
    }

    // Check OpenGL errors
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error before particle setup: " << err << std::endl;
    }

    // Generate and bind VAO
    glGenVertexArrays(1, &vertArrObj);
    glBindVertexArray(vertArrObj);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after generating VAO: " << err << std::endl;
        return;
    }

    // Generate and populate position buffer
    glGenBuffers(1, &vertBuffObj);
    glBindBuffer(GL_ARRAY_BUFFER, vertBuffObj);
    glBufferData(GL_ARRAY_BUFFER, sizeof(points), &points[0], GL_DYNAMIC_DRAW);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after setting up position buffer: " << err << std::endl;
    }

    // Generate and populate color buffer
    glGenBuffers(1, &colorBuffObj);
    glBindBuffer(GL_ARRAY_BUFFER, colorBuffObj);
    glBufferData(GL_ARRAY_BUFFER, sizeof(pointColors), &pointColors[0], GL_DYNAMIC_DRAW);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after setting up color buffer: " << err << std::endl;
    }

    // Unbind buffers
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    std::cout << "Particle system initialization complete" << std::endl;
    std::cout << "VAO: " << vertArrObj << ", Position VBO: " << vertBuffObj << ", Color VBO: " << colorBuffObj << std::endl;
}

void particleSys::reSet() {
	for (int i=0; i < numP; i++) {
		particles[i]->load(start);
	}
}

void particleSys::drawMe(std::shared_ptr<Program> prog) {
    // Check if program is valid
    if (!prog) {
        std::cerr << "ERROR: Null program passed to particleSys::drawMe" << std::endl;
        return;
    }

    // Check if vertex array object is initialized
    if (vertArrObj == 0) {
        std::cerr << "ERROR: Particle VAO not initialized. Did you call gpuSetup()?" << std::endl;
        return;
    }

    // Check OpenGL errors before we start
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error before particle drawing: " << err << std::endl;
    }

    // Try to bind the VAO and check for errors
    glBindVertexArray(vertArrObj);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after binding VAO: " << err << std::endl;
        return;
    }

    // Enable blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    // Try to enable vertex attributes
    glEnableVertexAttribArray(0); // For position
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after enabling position attribute: " << err << std::endl;
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Bind position buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertBuffObj);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after binding position buffer: " << err << std::endl;
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Set up position attribute pointer
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after setting position attribute pointer: " << err << std::endl;
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Set up instance divisor for position
    glVertexAttribDivisor(0, 1);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after setting position divisor: " << err << std::endl;
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Enable color attribute
    glEnableVertexAttribArray(1); // For color
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after enabling color attribute: " << err << std::endl;
        glVertexAttribDivisor(0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Bind color buffer
    glBindBuffer(GL_ARRAY_BUFFER, colorBuffObj);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after binding color buffer: " << err << std::endl;
        glDisableVertexAttribArray(1);
        glVertexAttribDivisor(0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Set up color attribute pointer
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after setting color attribute pointer: " << err << std::endl;
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(1);
        glVertexAttribDivisor(0, 0);
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Set up instance divisor for color
    glVertexAttribDivisor(1, 1);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after setting color divisor: " << err << std::endl;
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableVertexAttribArray(1);
        glVertexAttribDivisor(0, 0);
        glDisableVertexAttribArray(0);
        glBindVertexArray(0);
        glDisable(GL_BLEND);
        return;
    }

    // Attempt to draw
    glDrawArraysInstanced(GL_POINTS, 0, 1, numP);
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error after draw call: " << err << std::endl;
    }

    // Clean up
    glVertexAttribDivisor(1, 0);
    glVertexAttribDivisor(0, 0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(0);
    glBindVertexArray(0);
    glDisable(GL_BLEND);

    // Final error check
    err = glGetError();
    if (err != GL_NO_ERROR) {
        std::cerr << "OpenGL error at end of particle drawing: " << err << std::endl;
    }
}

void particleSys::setEmitter(const vec3& pos) {
	emitterPos = pos;
}

void particleSys::update() {

  vec3 pos;
  vec4 col;

  //update the particles
  for(auto particle : particles) {
        particle->update(t, h, g, emitterPos);
  }
  t += h;
 
  // Sort the particles by Z
  //temp->rotate(camRot, vec3(0, 1, 0));
  //be sure that camera matrix is updated prior to this update
  vec3 s, t, sk;
  vec4 p;
  quat r;
  glm::decompose(theCamera, s, r, t, sk, p);
  sorter.C = glm::toMat4(r); 
  sort(particles.begin(), particles.end(), sorter);


  //go through all the particles and update the CPU buffer
   for (int i = 0; i < numP; i++) {
        pos = particles[i]->getPosition();
        col = particles[i]->getColor();
		//update positions
        points[i*3+0] =pos.x; 
        points[i*3+1] =pos.y; 
        points[i*3+2] =pos.z; 
		//update colors
		pointColors[i * 4 + 0] = col.r;// +col.a / 10;
		pointColors[i * 4 + 1] = col.g;// +col.g / 10;
		pointColors[i * 4 + 2] = col.b;// +col.b / 10;
        pointColors[i * 4 + 3] = col.a;
  } 

  //update the GPU data
   glBindBuffer(GL_ARRAY_BUFFER, vertBuffObj);
   glBufferData(GL_ARRAY_BUFFER, sizeof(points), NULL, GL_STREAM_DRAW);
   glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*numP*3, points);
   //update gpu data for colors
   glBindBuffer(GL_ARRAY_BUFFER, colorBuffObj);
   glBufferData(GL_ARRAY_BUFFER, sizeof(pointColors), NULL, GL_STREAM_DRAW);
   glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*numP*4, pointColors);
   glBindBuffer(GL_ARRAY_BUFFER, 0);

}
