//
// sueda - geometry edits Z. Wood
// 3/16
//

#include <iostream>
#include "Particle.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Texture.h"


float randFloat(float l, float h)
{
	float r = rand() / (float) RAND_MAX;
	return (1.0f - r) * l + r * h;
}

Particle::Particle(vec3 start) :
	charge(1.0f),
	m(1.0f),
	d(0.0f),
	x(start),
	v(0.0f, 0.0f, 0.0f),
	lifespan(1.0f),
	tEnd(0.0f),
	scale(1.0f),
	color(1.0f, 1.0f, 1.0f, 1.0f)
{
}

Particle::~Particle()
{
}

void Particle::load(vec3 start)
{
	// Random initialization
	rebirth(0.0f, start);
}

/* all particles born at the origin */
void Particle::rebirth(float t, vec3 emitterPos)
{
	// sample a random point around the mesh center.
	float offsetX = randFloat(-0.5f, 0.5f);
	float offsetY = randFloat(-0.5f, 0.5f);
	vec3 randomOffset = vec3(offsetX, offsetY, 0.0);

	x = emitterPos + randomOffset;

	charge = randFloat(0.0f, 1.0f) < 0.5 ? -1.0f : 1.0f;	
	m = 1.0f;
  	d = randFloat(0.0f, 0.02f);
	//x = emitterPos;
	v.x = randFloat(-0.27f, 0.3f);
	v.y = randFloat(-0.1f, 0.9f);
	v.z = randFloat(-0.3f, 0.27f);
	lifespan = randFloat(1.0f, 5.0f); 
	tEnd = t + lifespan;
	scale = randFloat(0.2, 1.0f);
	// base color is yellow (RGB = 1,1,0) with full opacity
   	color.r = 1.0f + randFloat(-0.5f, 0.0f);
   	color.g = 1.0f + randFloat(-0.5f, 0.0f);
   	color.b = 0.0f + randFloat(-0.5f, 0.0f);
	color.a = 1.0f;
}

void Particle::update(float t, float h, const vec3 &g, const vec3 &emitterPos)
{
	if(t > tEnd) {
		rebirth(t, emitterPos);
	}

	// compute a radial direction from the emitter to the particle
	vec3 radialVec = x - emitterPos;
	// avoid division by zero if the particle is at the emitter
	vec3 radialDir = (length(radialVec) > 0.001f) ? normalize(radialVec) : vec3(0.0f);

	// radial acceleration: push the particle away from the emitter
	float radialAccelFactor = 5.0f;  // creates more explosive acceleration
	vec3 radialAccel = radialAccelFactor * radialDir;

	// prof said random = real life so add a bit of random turbulence for a natural look
	vec3 randomAccel(randFloat(-0.05f, 0.05f),
		randFloat(-0.05f, 0.05f),
		randFloat(-0.05f, 0.05f));

	// total acceleration: gravity + radial force + turbulence
	vec3 acceleration = radialAccel + randomAccel;

	// update the velocity with the acceleration
	v += acceleration * h;

	// apply damping to simulate drag (prevent runaway speeds!!)
	float damping = 0.99f; // value slightly less than 1.0
	v *= damping;

	// update position based on the new velocity.
	x += v * h;

	//very simple update
	//x += h*v;
	//To do - how do you want to update the forces?

	// update color opacity based on lifespan remaining
	color.a = (tEnd-t)/lifespan;
}
