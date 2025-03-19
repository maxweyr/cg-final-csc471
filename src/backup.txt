/*
 * Program 4 base code - includes modifications to shape and initGeom in preparation to load
 * multi shape objects 
 * CPE 471 Cal Poly Z. Wood + S. Sueda + I. Dunn
 * Max Eyrich
 */

#include <iostream>
#include <glad/glad.h>
#include "GLSL.h"
#include "Program.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "WindowManager.h"
#include "Texture.h"
#define _USE_MATH_DEFINES
#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader/tiny_obj_loader.h>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <unordered_map>
#include <random>
#include "stb_image.h"
#include "Bezier.h"
#include "Spline.h"
#include "particleSys.h"

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

using namespace std;
using namespace glm;

class Application : public EventCallbacks {
public:

	WindowManager * windowManager = nullptr;

	particleSys* doorParticleSystem;

	// shader programs: without texture, with textures
	std::shared_ptr<Program> prog, texProg, partProg;

	// geometry
	std::vector<std::shared_ptr<Shape>> meshes;

	struct MultiMesh {
		vector<shared_ptr<Shape>> shapes;  // Each sub-mesh in the object
		float scale;                       // Uniform scale for this object
	};

	unordered_map<string, MultiMesh> multiMeshes;  // Store multi-mesh objects

	// containers for static placements
	std::vector<glm::mat4> starPlacements;
	unordered_map<string, vector<glm::mat4>> rockPlacements;
	unordered_map<string, vector<vector<glm::mat4>>> foliagePlacements, treePlacements;

	// global data for ground plane - direct load constant defined CPU data to GPU (not obj)
	GLuint GrndBuffObj, GrndNorBuffObj, GrndTexBuffObj, GIndxBuffObj;
	int g_GiboLen;
	// ground VAO
	GLuint GroundVertexArrayID;

	// images to be used as textures
	shared_ptr<Texture> texture0, texture1, texture2, particleTexture; // ground, sky, door, particle textures

	// global data
	vec3 gMin;
	float gRot = 0;
	float gCamH = 0;

	// camera variables
	float phi = 0.0f;								// Pitch angle in radians
	float theta = glm::pi<float>();					// Yaw angle in radians (start facing backward to match my current view)
	float radius = 15.0f;							// Distance from camera to look-at point
	glm::vec3 eye = glm::vec3(0.0f, 0.0f, 0.0f);	// Camera position
	glm::vec3 lookAt;								// Point to look at
	double lastX, lastY;
	bool mouseInitialized = false;
	float mouseSpeed = 0.03f;						// sensitivty

	// For cinematic camera tour
	bool tourEnabled = false;              // Flag to indicate tour mode
	Spline* cameraTour = nullptr;          // Pointer to a spline for the tour
	glm::vec3 tourLookAt = glm::vec3(0.0f, 0.0f, 0.0f); // Fixed look-at point for the tour
	float tourDuration = 10.0f;            // Duration (in seconds) for the tour
	float lastTourUpdateTime = 0.0f;       // For computing delta time in tour mode

	float groundSize = 20.0f;
	float sphereRadius = groundSize / 2.0f;

	// star perameters
	int numStars = 100;						// number of stars
	std::vector<float> starBaseScales;		// store original star scales
	std::vector<glm::vec3> starPositions;	// store star positions
	std::vector<float> starPhases;			// random phase offset for each star
	std::vector<float> starFrequencies;		// different frequency for each star
	std::vector<float> starAlphas;			// current alpha/opacity of each star
	std::vector<float> starFadePhases;		// different phase for fade animation
	std::vector<float> starFadeRates;		// how quickly stars fade in/out

	int treeIndices[3] = { 11, 13, 15 };			// specific meshes from a multi mesh

	// animation data
	float moonSpeed = 0.15f;						// angular speed (rad/s)
	float moonY = sphereRadius * 0.66f;				// orbit height

	// particle parameters
	bool keyToggles[256] = { false };
	float t = 0.0f; //reset in init
	float h = 0.01f;
	glm::vec3 doorPos = glm::vec3(0.0f, 0.0f, ((sphereRadius * multiMeshes["sphereWTex"].scale * 2.0f) - 1.0f)); // place door at behind starting camera orientation

	void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
	{
		if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) { // exit the program when 'ESC' is pressed
			// shutdown particle system
			delete doorParticleSystem;
			glfwSetWindowShouldClose(window, GL_TRUE);
		}

		if (key == GLFW_KEY_G && action == GLFW_PRESS) { // toggle cinematic tour when 'g' is pressed
			if (!tourEnabled) { // start tour
				tourLookAt = glm::vec3(0.0f, 0.0f, 0.0f); // fix the look-at point
				glm::vec3 start = vec3(15, 5, 0); // use current eye position as the start of the tour
				glm::vec3 control = start +
					glm::vec3(0.0f, 5.0f, 15.0f); // control point to create a nice arc
				glm::vec3 end = glm::vec3(-15.0f, 5.0f, 0.0f); // end point (destination) for camera to travel to
				cameraTour = new Spline(start, control,
					end, tourDuration); // create a quadratic spline (order2) that takes tourDuration seconds
				tourEnabled = true;
				lastTourUpdateTime = glfwGetTime(); // reset timer for the spline update
			} else {
				tourEnabled = false; // turn off tour mode
				if (cameraTour) {
					delete cameraTour;
					cameraTour = nullptr;
				}
				theta = glm::pi<float>(); // reset camera parameters for manual control
				phi = 0.0f;
				radius = 15.0f;
				updateLookAt();
			}
		}

		if (!tourEnabled) {	// if not in tour mode allow manual movement
			if (action == GLFW_PRESS || action == GLFW_REPEAT) {
				glm::vec3 forward = glm::normalize(lookAt - eye);
				glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
				if (key == GLFW_KEY_W)
					eye += forward * mouseSpeed;
				if (key == GLFW_KEY_S)
					eye -= forward * mouseSpeed;
				if (key == GLFW_KEY_A)
					eye -= right * mouseSpeed;
				if (key == GLFW_KEY_D)
					eye += right * mouseSpeed;
				updateLookAt();
			}
		}
		if (key == GLFW_KEY_Z && action == GLFW_PRESS) { // geometry view when 'z' pressed
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		}
		if (key == GLFW_KEY_Z && action == GLFW_RELEASE) { // undo geometry view when 'z' released
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		}
	}

	void mouseMoveCallback(GLFWwindow* window, double xpos, double ypos) {
		if (!mouseInitialized) {
			lastX = xpos;
			lastY = ypos;
			mouseInitialized = true;
			return;
		}

		// Calculate the change in mouse position
		float deltaX = xpos - lastX;
		float deltaY = ypos - lastY;
		lastX = xpos;
		lastY = ypos;

		// update theta (yaw) based on x movement
		theta += deltaX * mouseSpeed;
		// keep theta in [0, 2π]
		if (theta > 2 * glm::pi<float>()) theta -= 2 * glm::pi<float>();
		if (theta < 0) theta += 2 * glm::pi<float>();

		// update phi (pitch) based on y movement
		phi -= deltaY * mouseSpeed;

		// constrain phi to prevent looking too far up or down (80 degrees up or down)
		const float maxPhi = 80.0f * glm::pi<float>() / 180.0f;  // 80 degrees in radians
		phi = glm::clamp(phi, -maxPhi, maxPhi);

		// calculate the new lookAt point based on phi and theta
		updateLookAt();
	}

	void scrollCallback(GLFWwindow* window, double deltaX, double deltaY) {
		radius -= deltaY * 0.5f;						// adjust radius based on scroll (zoom in/out)
		radius = glm::clamp(radius, 1.0f, 50.0f);		// keep radius within reasonable bounds
		updateLookAt();									// update lookAt point with new radius
	}

	void updateLookAt() {
		// calculate the lookAt point using spherical to Cartesian conversion
		lookAt.x = eye.x + radius * cos(phi) * cos(theta);
		lookAt.y = eye.y + radius * sin(phi);
		lookAt.z = eye.z + radius * cos(phi) * cos((glm::pi<float>() / 2.0) - theta);
	}

	void mouseCallback(GLFWwindow* window, int button, int action, int mods) override {
		std::cout << "Mouse button " << button << " action " << action << std::endl;
	}

	void resizeCallback(GLFWwindow *window, int width, int height) {
		glViewport(0, 0, width, height);
	}

	void init(const std::string& resourceDirectory) {
		GLSL::checkVersion();
		glClearColor(.72f, .84f, 1.06f, 1.0f);		// set background color
		glEnable(GL_DEPTH_TEST);					// enable z-buffer test

		// Initialize the GLSL program that we will use for local shading
		prog = make_shared<Program>();
		prog->setVerbose(true);
		prog->setShaderNames(resourceDirectory + "/simple_vert.glsl", resourceDirectory + "/simple_frag.glsl");
		prog->init();
		prog->addUniform("P");
		prog->addUniform("V");
		prog->addUniform("M");
		prog->addUniform("flip");
		prog->addUniform("MatAmb");
		prog->addUniform("MatDif");
		prog->addUniform("MatSpec");
		prog->addUniform("MatShine");
		prog->addUniform("lightPos");
		prog->addUniform("lightIntensity");
		prog->addUniform("uTime");
		prog->addUniform("windStrength");
		prog->addUniform("windFrequency");
		prog->addUniform("applyWind");
		prog->addUniform("ambientLight");
		prog->addUniform("useFixedLight");
		prog->addUniform("fixedLightPos");
		prog->addUniform("ignoreLight");
		prog->addUniform("starColor");
		prog->addAttribute("vertPos");
		prog->addAttribute("vertNor");
		prog->addAttribute("vertTex");

		// Initialize the GLSL program that we will use for texture mapping
		texProg = make_shared<Program>();
		texProg->setVerbose(true);
		texProg->setShaderNames(resourceDirectory + "/tex_vert.glsl", resourceDirectory + "/tex_frag.glsl");
		texProg->init();
		texProg->addUniform("P");
		texProg->addUniform("V");
		texProg->addUniform("M");
		texProg->addUniform("flip");
		texProg->addUniform("Texture");
		texProg->addUniform("MatShine");
		texProg->addUniform("lightPos");
		texProg->addAttribute("vertPos");
		texProg->addAttribute("vertNor");
		texProg->addAttribute("vertTex");

		partProg = make_shared<Program>();
		partProg->setVerbose(true);
		partProg->setShaderNames(resourceDirectory + "/part_vert.glsl", resourceDirectory + "/part_frag.glsl");
		if (!partProg->init()) {
			std::cerr << "Particle shader failed to compile... exiting!" << std::endl;
			exit(1);
		}
		partProg->addUniform("P");
		partProg->addUniform("M");
		partProg->addUniform("V");
		partProg->addUniform("pColor");
		partProg->addUniform("alphaTexture");
		partProg->addAttribute("vertPos");

		// create particle system (set the emitter position later)
		doorParticleSystem = new particleSys(vec3(0, 0, 0));
		doorParticleSystem->gpuSetup();
	}

	void initGeom(const std::string& resourceDirectory) {
		// mesh objects to target
		vector<string> multiMeshFiles = {
			"castleDoor.obj",
			"foliageCollection.obj",
			"moon.obj",
			"natureCollection.obj",
			"pottedFern.obj",
			"rock1.obj",
			"rock2.obj",
			"rock3.obj",
			"rock4.obj",
			"rock5.obj",
			"sphereWTex.obj",
			"star.obj",
			"cube.obj",
			"spaceKit.obj",
			"moonDoor.obj"
		};

		for (const auto& filename : multiMeshFiles) {
			vector<tinyobj::shape_t> shapes;
			vector<tinyobj::material_t> objMaterials;
			string errStr;

			bool rc = tinyobj::LoadObj(shapes, objMaterials, errStr, (resourceDirectory + "/models/" + filename).c_str());
			if (!rc) {
				cerr << "Failed to load " << filename << ": " << errStr << endl;
				continue;
			}

			MultiMesh multiMesh;
			glm::vec3 gMin(FLT_MAX), gMax(-FLT_MAX);

			// process each shape in the obj file
			for (auto& shape : shapes) {
				auto mesh = make_shared<Shape>();
				mesh->createShape(shape);
				mesh->measure();
				mesh->init();
				multiMesh.shapes.push_back(mesh);
				// update bounding box for scaling
				gMin = glm::min(gMin, mesh->min);
				gMax = glm::max(gMax, mesh->max);
			}

			// compute uniform scale factor
			float maxExtent = glm::max(gMax.x - gMin.x, glm::max(gMax.y - gMin.y, gMax.z - gMin.z));
			multiMesh.scale = 2.0f / maxExtent; // scale objects to fit within a reasonable size

			// use filesystem to get the stem
			fs::path filePath(filename);
			string key = filePath.stem().string();

			// store in the map using the simplified key
			multiMeshes[key] = multiMesh;
		}
	}

	void initStars() {
		// make sure the sphere mesh is loaded
		std::string sphereKey = "sphereWTex";
		if (multiMeshes.find(sphereKey) == multiMeshes.end() || multiMeshes[sphereKey].shapes.empty()) {
			std::cerr << "Missing sphere for star placement!" << std::endl;
			return;
		}

		// make sure the star mesh is loaded
		std::string starKey = "star";
		if (multiMeshes.find(starKey) == multiMeshes.end() || multiMeshes[starKey].shapes.empty()) {
			std::cerr << "Missing star mesh for star placement!" << std::endl;
			return;
		}

		starPlacements.resize(numStars); // use numStars to adjust star count
		std::default_random_engine rng(std::random_device{}());

		float maxHeight = sphereRadius * 0.85f; // ceiling height (85% sphere radius)
		float minDistance = sphereRadius * 0.12f; // mini distance between stars (12% sphere radius)

		// generate positions using Fibonacci sphere algorithm for better distribution
		const float goldenRatio = (1.0f + std::sqrt(5.0f)) / 2.0f;
		const float angleIncrement = 2.0f * glm::pi<float>() * goldenRatio;

		// keep track of placed stars
		std::vector<glm::vec3> placedPositions;
		std::vector<glm::mat4> validPlacements;
		std::vector<float> validScales;      // Store scales for animation
		std::vector<glm::vec3> validPositions; // Store positions for animation

		// generate more positions than needed for filtering
		int attempts = numStars * 3;
		std::uniform_real_distribution<float> scaleDist(0.5f, 1.0f); // random star scale
		std::uniform_real_distribution<float> randomOffsetDist(-0.05f, 0.05f); // small random offset

		for (int i = 0; i < attempts && validPlacements.size() < numStars; i++) {
			// using Fibonacci sphere for initial position
			float t = static_cast<float>(i) / attempts;
			float inclination = std::acos(1.0f - 2.0f * t);
			float azimuth = angleIncrement * i;

			// convert to Cartesian coordinates
			float R = sphereRadius * multiMeshes[sphereKey].scale * 0.95f; // small inner radius offset
			float x = R * std::sin(inclination) * std::cos(azimuth);
			float y = R * std::cos(inclination);
			float z = R * std::sin(inclination) * std::sin(azimuth);

			// only place stars in upper hemisphere with respect to maxHeight
			if (y > 0 && y <= maxHeight) {
				// add small random offset to break up the perfect pattern
				x += R * randomOffsetDist(rng);
				z += R * randomOffsetDist(rng);

				glm::vec3 starPos(x, y, z);

				// check distance from all previously placed stars
				bool tooClose = false;
				for (const auto& pos : placedPositions) {
					if (glm::distance(pos, starPos) < minDistance) {
						tooClose = true;
						break;
					}
				}

				if (!tooClose) {
					// compute rotation so star is facing the origin
					glm::mat4 rotation = glm::inverse(glm::lookAt(
						starPos,
						glm::vec3(0.0f),
						glm::vec3(0.0f, 1.0f, 0.0f)
					));

					glm::mat4 translation = glm::translate(glm::mat4(1.0f), starPos);
					float scale = multiMeshes[starKey].scale * scaleDist(rng);
					glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(scale));

					// combine transforms
					glm::mat4 model = translation * rotation * scaling;

					validPlacements.push_back(model);
					placedPositions.push_back(starPos);

					// store values needed for animation
					validScales.push_back(scale);
					validPositions.push_back(starPos);
				}
			}
		}

		// if all star placements aren't filled, fill remaining with random positions
		if (validPlacements.size() < numStars) {
			std::uniform_real_distribution<float> heightDist(0.0f, maxHeight);
			std::uniform_real_distribution<float> angleDist(0.0f, glm::two_pi<float>());

			while (validPlacements.size() < numStars) {
				float theta = angleDist(rng);
				float phi = angleDist(rng) / 2.0f; // only use top half of the sphere
				float h = heightDist(rng);
				float R = sphereRadius * multiMeshes[sphereKey].scale * 0.95f;
				float rSlice = sqrt(R * R - h * h);

				float x = rSlice * cos(theta);
				float z = rSlice * sin(theta);
				glm::vec3 starPos(x, h, z);

				glm::mat4 rotation = glm::inverse(glm::lookAt(
					starPos,
					glm::vec3(0.0f),
					glm::vec3(0.0f, 1.0f, 0.0f)
				));

				glm::mat4 translation = glm::translate(glm::mat4(1.0f), starPos);
				float scale = multiMeshes[starKey].scale * scaleDist(rng);
				glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(scale));

				glm::mat4 model = translation * rotation * scaling;
				validPlacements.push_back(model);

				// store values needed for animation
				validScales.push_back(scale);
				validPositions.push_back(starPos);
			}
		}

		// copy valid placements into starPlacements vector
		std::copy(validPlacements.begin(), validPlacements.end(), starPlacements.begin());

		// initialize animation data
		starBaseScales.resize(numStars);
		starPositions.resize(numStars);
		starPhases.resize(numStars);
		starFrequencies.resize(numStars);

		std::uniform_real_distribution<float> phaseDist(0.0f, glm::two_pi<float>());
		std::uniform_real_distribution<float> freqDist(0.5f, 1.5f);

		for (int i = 0; i < numStars; i++) {
			starBaseScales[i] = validScales[i];
			starPositions[i] = validPositions[i];
			starPhases[i] = phaseDist(rng);
			starFrequencies[i] = freqDist(rng);
		}

		initStarFade();

		std::cout << "Successfully placed " << starPlacements.size() << " stars" << std::endl;
	}

	void initStarFade() {
		// Initialize fade parameters
		starAlphas.resize(numStars, 1.0f);        // Start fully visible
		starFadePhases.resize(numStars);
		starFadeRates.resize(numStars);

		std::default_random_engine rng(std::random_device{}());
		std::uniform_real_distribution<float> phaseDist(0.0f, glm::two_pi<float>());
		std::uniform_real_distribution<float> rateDist(0.2f, 0.6f);   // Slower fade rate

		for (int i = 0; i < numStars; i++) {
			starFadePhases[i] = phaseDist(rng);
			starFadeRates[i] = rateDist(rng);
		}
	}

	void updateStarAnimation(float currentTime) {
		for (int i = 0; i < numStars; i++) {
			// calculate new scale with pulsing effect
			float pulseAmount = 0.15f; // 15% size variation
			float pulseFactor = 1.0f + pulseAmount * sin(starFrequencies[i] * currentTime + starPhases[i]);
			float newScale = starBaseScales[i] * pulseFactor;
			// calculate new fade to enhance pulsing
			float fadeFreq = starFadeRates[i];
			float fadeVal = sin(fadeFreq * currentTime + starFadePhases[i]);

			// map from [-1, 1] to [0.0, 1.0] for alpha
			starAlphas[i] = (fadeVal * 0.5f + 0.5f);

			// new matrix using the original position and orientation
			glm::vec3 starPos = starPositions[i];

			// calculate rotation the same way as in initStars
			glm::mat4 rotation = glm::inverse(glm::lookAt(
				starPos,
				glm::vec3(0.0f),
				glm::vec3(0.0f, 1.0f, 0.0f)
			));

			glm::mat4 translation = glm::translate(glm::mat4(1.0f), starPos);
			glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(newScale));

			// update the star's placement
			starPlacements[i] = translation * rotation * scaling;
		}
	}

	void initFoliage() {
		// Use the key for your foliage collection (make sure it matches the filename used in initGeom)
		string foliageKey = "foliageCollection";
		if (multiMeshes.find(foliageKey) == multiMeshes.end()) {
			cerr << "Foliage collection not loaded." << endl;
			return;
		}
		MultiMesh& foliage = multiMeshes[foliageKey];
		vector<vector<glm::mat4>> placements;	// prepare a container for the placements
		placements.resize(foliage.shapes.size());	// one vector per plant type

		// define inner and outer radii for annular region
		float rInner = sphereRadius * multiMeshes["sphereWTex"].scale / 2.0f;
		float rOuter = sphereRadius * multiMeshes["sphereWTex"].scale * 2.0f;

		// setup random number generators
		std::default_random_engine rng(std::random_device{}());
		std::uniform_int_distribution<int> countDist(50, 100);						// random range of plants per type
		std::uniform_real_distribution<float> thetaDist(0.0f, glm::two_pi<float>());
		std::uniform_real_distribution<float> rDist(rInner, rOuter);				// random radius
		std::uniform_real_distribution<float> rotDist(0.0f, glm::two_pi<float>());	// random rotation about Y

		// For each plant (sub-mesh) in the foliage collection:
		for (size_t i = 0; i < foliage.shapes.size(); i++) {
			int instanceCount = countDist(rng);
			placements[i].resize(instanceCount);
			for (int j = 0; j < instanceCount; j++) {
				glm::mat4 model = glm::mat4(1.0f);									// build the model matrix:
				float r = rDist(rng);
				float theta = thetaDist(rng);
				float x = r * cos(theta);
				float z = r * sin(theta);
				float y = -1.25f;													// fix the vertical position
				float angle = rotDist(rng);											// random rotation about Y
				model = glm::translate(model, glm::vec3(x, y, z));					// set translation
				model = glm::rotate(model, angle, glm::vec3(0, 1, 0));				// set random rotation
				model = glm::scale(model, glm::vec3(foliage.scale));				// set precomputed uniform scale factor
				placements[i][j] = model;
			}
		}
		foliagePlacements[foliageKey] = placements;
	}

	void Application::initNature() {
		std::string natureKey = "natureCollection";
		if (multiMeshes.find(natureKey) == multiMeshes.end()) {
			std::cerr << "Nature collection not loaded." << std::endl;
			return;
		}
		MultiMesh& nature = multiMeshes[natureKey];
		vector<vector<glm::mat4>> placements; // prepare a container for the placements
		placements.resize(3); // one vector per tree type

		// define inner and outer radii for annular region
		float rInner = 4.0f * sphereRadius * multiMeshes["sphereWTex"].scale / 3.0f;
		float rOuter = 2.0f * sphereRadius * multiMeshes["sphereWTex"].scale - 5.0f;

		// set up random generators
		std::default_random_engine rng(std::random_device{}());
		std::uniform_real_distribution<float> rDist(rInner, rOuter);
		std::uniform_real_distribution<float> thetaDist(0.0f, glm::two_pi<float>());
		std::uniform_int_distribution<int> countDistT(50, 70); // random number of trees
		std::uniform_int_distribution<int> countDistLS(10, 20); // random number of stumps and logs
		std::uniform_real_distribution<float> rotDist(0.0f, glm::two_pi<float>()); // random rotation

		// define door position and the path that should be kept clear
		float sphereEdge = sphereRadius * multiMeshes["sphereWTex"].scale * 2.0f;
		glm::vec3 doorPos(0, -1.25f, sphereEdge - 1.0f); // door position
		glm::vec3 centerPos(0, -1.25f, 0); // center of the scene

		// calculate door-to-center direction vector
		glm::vec3 doorToCenter = glm::normalize(centerPos - doorPos);

		// define a corridor width (adjust as needed)
		float corridorWidth = 5.0f;  // width of clear path

		for (int i = 0; i < 3; i++) { // range of placements per type
			int instanceCount = i < 1 ? countDistT(rng) : countDistLS(rng);
			// vector to store valid placements for this tree type
			std::vector<glm::mat4> validPlacements;

			// try to place more trees than needed, then filter out those in the corridor
			int attemptsLimit = instanceCount * 2;  // generate twice as many candidates

			for (int j = 0; j < attemptsLimit && validPlacements.size() < instanceCount; j++) {
				// choose a random radius and angle within the annular region
				float r = rDist(rng);
				float theta = thetaDist(rng);
				float x = r * cos(theta);
				float z = r * sin(theta);
				float y = -1.25f; // fix the vertical position

				// create potential position
				glm::vec3 treePos(x, y, z);

				// check if the tree would block the path from door to center
				// calculate distance from point to line (door-to-center)
				glm::vec3 doorToTree = treePos - doorPos;
				float projectionLength = glm::dot(doorToTree, doorToCenter);

				// only consider trees that are between the door and center (not behind the door)
				bool inPathDirection = projectionLength > 0 && projectionLength < glm::length(centerPos - doorPos);

				// calculate perpendicular distance from tree to the door-center line
				glm::vec3 projection = doorPos + doorToCenter * projectionLength;
				float distanceToPath = glm::length(treePos - projection);

				// skip this position if the tree would be in the corridor
				if (inPathDirection && distanceToPath < corridorWidth) {
					continue;  // Skip this position
				}

				// position is valid, create the model matrix
				float angle = rotDist(rng); // random rotation about Y
				glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z)); // build model matrix and set translation
				model = glm::rotate(model, angle, glm::vec3(0, 1, 0)); // set random rotation
				model = glm::scale(model, 5.0f * glm::vec3(1, 1, 1)); // set scale

				validPlacements.push_back(model);
			}

			// store the valid placements (might be less than requested if many were filtered out)
			placements[i] = validPlacements;
		}

		treePlacements[natureKey] = placements;
	}

	void initRocks() {
		vector<string> rockKeys = { "rock1", "rock2", "rock3", "rock4", "rock5" };

		// define inner and outer radii for annular region
		float rInner = 4.0f * sphereRadius * multiMeshes["sphereWTex"].scale / 3.0f;
		float rOuter = 2.0f * sphereRadius * multiMeshes["sphereWTex"].scale;

		// Set up random generators.
		std::default_random_engine rng(std::random_device{}());
		std::uniform_int_distribution<int> countDist(10, 30);  // random number of instances per rock type
		std::uniform_real_distribution<float> thetaDist(0.0f, glm::two_pi<float>());
		std::uniform_real_distribution<float> rDist(rInner, rOuter);
		std::uniform_real_distribution<float> rotDist(0.0f, glm::two_pi<float>());
		std::uniform_real_distribution<float> scaleDist(1.0f, 1.5f); // random scale variation for rocks

		for (auto& rockKey : rockKeys) {	// for each rock type, generate a set of placements
			if (multiMeshes.find(rockKey) == multiMeshes.end())	// check existence
				continue;

			auto rockShape = multiMeshes[rockKey].shapes[0]; // get the rock's shape.
			float rockYOffset = -1.5f - rockShape->min.y;	// compute an offset so the rock's lowest point touches the ground

			int instanceCount = countDist(rng);
			vector<glm::mat4> placements;
			placements.resize(instanceCount);

			for (int i = 0; i < instanceCount; i++) {	// for each rock choose a random radius and angle within the annular region
				float r = rDist(rng);
				float theta = thetaDist(rng);
				float x = r * cos(theta);
				float z = r * sin(theta);
				glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(x, rockYOffset, z));	// build the model matrix and translate
				model = glm::rotate(model, rotDist(rng), glm::vec3(0, 1, 0));	// apply a random rotation about the Y axis
				model = glm::scale(model, glm::vec3(scaleDist(rng)));			// apply random scaling
				placements[i] = model;
			}
			rockPlacements[rockKey] = placements;
		}
	}

	// code to load in the new textures
	void initTex(const std::string& resourceDirectory) {
		texture0 = make_shared<Texture>();
		texture0->setFilename(resourceDirectory + "/cartoonGrass.jpg");
		texture0->init();
		texture0->setUnit(0);
		texture0->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);
		texture1 = make_shared<Texture>();
		texture1->setFilename(resourceDirectory + "/starSky.jpg");
		texture1->init();
		texture1->setUnit(1);
		texture1->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);
		texture2 = make_shared<Texture>();
		texture2->setFilename(resourceDirectory + "/bricks.jpg");
		texture2->init();
		texture2->setUnit(2);
		texture2->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);
		particleTexture = make_shared<Texture>();
		particleTexture->setFilename(resourceDirectory + "/alpha.bmp");
		particleTexture->init();
		particleTexture->setUnit(3);
		particleTexture->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);
	}

	//directly pass quad for the ground to the GPU
	void initGround() {
		float g_groundSize = groundSize;
		float g_groundY = -0.25;
  		// A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
		float GrndPos[] = {
			-g_groundSize, g_groundY, -g_groundSize,
			-g_groundSize, g_groundY,  g_groundSize,
			g_groundSize, g_groundY,  g_groundSize,
			g_groundSize, g_groundY, -g_groundSize
		};
		float GrndNorm[] = {
			0, 1, 0,
			0, 1, 0,
			0, 1, 0,
			0, 1, 0,
			0, 1, 0,
			0, 1, 0
		};
		static GLfloat GrndTex[] = {
      		0, 0, // back
      		0, 1,
      		1, 1,
      		1, 0
		};
      	unsigned short idx[] = {0, 1, 2, 0, 2, 3};
		//generate the ground VAO
      	glGenVertexArrays(1, &GroundVertexArrayID);
      	glBindVertexArray(GroundVertexArrayID);

      	g_GiboLen = 6;
      	glGenBuffers(1, &GrndBuffObj);
      	glBindBuffer(GL_ARRAY_BUFFER, GrndBuffObj);
      	glBufferData(GL_ARRAY_BUFFER, sizeof(GrndPos), GrndPos, GL_STATIC_DRAW);

      	glGenBuffers(1, &GrndNorBuffObj);
      	glBindBuffer(GL_ARRAY_BUFFER, GrndNorBuffObj);
      	glBufferData(GL_ARRAY_BUFFER, sizeof(GrndNorm), GrndNorm, GL_STATIC_DRAW);

      	glGenBuffers(1, &GrndTexBuffObj);
      	glBindBuffer(GL_ARRAY_BUFFER, GrndTexBuffObj);
      	glBufferData(GL_ARRAY_BUFFER, sizeof(GrndTex), GrndTex, GL_STATIC_DRAW);

      	glGenBuffers(1, &GIndxBuffObj);
     	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GIndxBuffObj);
      	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(idx), idx, GL_STATIC_DRAW);
      }

      //code to draw the ground plane
     void drawGround(shared_ptr<Program> curS) {
     	curS->bind();
     	glBindVertexArray(GroundVertexArrayID);
     	texture0->bind(curS->getUniform("Texture"));
		//draw the ground plane 
  		SetModel(vec3(0, -1, 0), 0, 0, 1, curS);
  		glEnableVertexAttribArray(0);
  		glBindBuffer(GL_ARRAY_BUFFER, GrndBuffObj);
  		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

  		glEnableVertexAttribArray(1);
  		glBindBuffer(GL_ARRAY_BUFFER, GrndNorBuffObj);
  		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

  		glEnableVertexAttribArray(2);
  		glBindBuffer(GL_ARRAY_BUFFER, GrndTexBuffObj);
  		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0);

   		// draw!
  		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GIndxBuffObj);
  		glDrawElements(GL_TRIANGLES, g_GiboLen, GL_UNSIGNED_SHORT, 0);

  		glDisableVertexAttribArray(0);
  		glDisableVertexAttribArray(1);
  		glDisableVertexAttribArray(2);
  		curS->unbind();
     }

     //helper function to pass material data to the GPU
	void SetMaterial(shared_ptr<Program> curS, int i) {

    	switch (i) {
    		case 0: // flat grey
				glUniform3f(curS->getUniform("MatAmb"), 0.2, 0.2, 0.2);
				glUniform3f(curS->getUniform("MatDif"), 0.5, 0.5, 0.5);
				glUniform3f(curS->getUniform("MatSpec"), 0.3, 0.3, 0.3);
				glUniform1f(curS->getUniform("MatShine"), 10.0);
    			break;
			case 1: // cartoonish moon material provided by chatGPT
				glUniform3f(curS->getUniform("MatAmb"), 0.1, 0.1, 0.1);		// high ambient for a luminous effect
				glUniform3f(curS->getUniform("MatDif"), 1.0f, 0.97f, 0.9f); // diffuse: matching off-white color
				glUniform3f(curS->getUniform("MatSpec"), 0.2f, 0.2f, 0.15f);// specular: low for a flat, cartoon feel
				glUniform1f(curS->getUniform("MatShine"), 8.0f);            // lower shininess for soft highlights
				break;
    		case 2: // cartoonish yellowish green material for grass and plants by chatGPT
				glUniform3f(curS->getUniform("MatAmb"), 0.05f, 0.06f, 0.035f);  // ambient: muted yellow-green
				glUniform3f(curS->getUniform("MatDif"), 0.2f, 0.6f, 0.2f);   // diffuse: vibrant, leafy green
				glUniform3f(curS->getUniform("MatSpec"), 0.1f, 0.1f, 0.1f);  // specular: very low to keep it flat
				glUniform1f(curS->getUniform("MatShine"), 5.0f);             // low shininess for a soft, cartoon feel
				break;
			case 3: // cartoonish yellow star material by chatGPT
				glUniform3f(curS->getUniform("MatAmb"), 0.25f, 0.20f, 0.05f);   // ambient: stronger yellow tint with less blue
				glUniform3f(curS->getUniform("MatDif"), 1.0f, 0.9f, 0.4f);     // diffuse: more saturated yellow with less blue/green
				glUniform3f(curS->getUniform("MatSpec"), 0.4f, 0.3f, 0.1f);    // specular: yellowish highlights
				glUniform1f(curS->getUniform("MatShine"), 20.0f);             // keep the same shininess
				break;
			case 4: // cartoonish purple star material by chatGPT
				glUniform3f(curS->getUniform("MatAmb"), 0.05f, 0.02f, 0.07f);   // ambient: moderate purple
				glUniform3f(curS->getUniform("MatDif"), 0.6f, 0.3f, 0.8f);   // diffuse: a bit brighter, still purple
				glUniform3f(curS->getUniform("MatSpec"), 0.2f, 0.1f, 0.3f);  // specular: low for a soft, cartoonish look
				glUniform1f(curS->getUniform("MatShine"), 15.0f);            // moderate shininess for subtle highlights
				break;
			case 5: // wood material for stumps, trunks, and logs
				glUniform3f(curS->getUniform("MatAmb"), 0.02f, 0.01f, 0.005f);  // ambient: dark, warm brown
				glUniform3f(curS->getUniform("MatDif"), 0.6f, 0.3f, 0.1f);   // diffuse: rich brown tone
				glUniform3f(curS->getUniform("MatSpec"), 0.2f, 0.2f, 0.2f);  // specular: modest highlights
				glUniform1f(curS->getUniform("MatShine"), 10.0f);            // moderate shininess
				break;
			case 6: // green material for leaves
				glUniform3f(curS->getUniform("MatAmb"), 0.01f, 0.03f, 0.01f);   // ambient: deep green
				glUniform3f(curS->getUniform("MatDif"), 0.2f, 0.6f, 0.2f);   // diffuse: vibrant, leafy green
				glUniform3f(curS->getUniform("MatSpec"), 0.1f, 0.1f, 0.1f);  // specular: very low to keep it flat
				glUniform1f(curS->getUniform("MatShine"), 5.0f);             // low shininess for a soft, cartoon feel
				break;
			case 7: // wanna be sandstone door material
				glUniform3f(curS->getUniform("MatAmb"), 0.08f, 0.07f, 0.06f);  // ambient: warm, light brown
				glUniform3f(curS->getUniform("MatDif"), 0.45f, 0.4f, 0.35f);  // diffuse: slightly brighter sandstone tone
				glUniform3f(curS->getUniform("MatSpec"), 0.2f, 0.2f, 0.2f);  // specular: low for subtle shine
				glUniform1f(curS->getUniform("MatShine"), 12.0f);           // moderate shininess for a soft highlight
				break;
			case 8: // shadow material: dark, flat
				glUniform3f(curS->getUniform("MatAmb"), 0.0f, 0.0f, 0.0f);
				glUniform3f(curS->getUniform("MatDif"), 0.0f, 0.0f, 0.0f);
				glUniform3f(curS->getUniform("MatSpec"), 0.0f, 0.0f, 0.0f);
				glUniform1f(curS->getUniform("MatShine"), 1.0f);
				break;
			case 9: // stone material for rocks
				glUniform3f(curS->getUniform("MatAmb"), 0.025f, 0.025f, 0.025f);  // ambient: dark grey
				glUniform3f(curS->getUniform("MatDif"), 0.6f, 0.6f, 0.6f);    // diffuse: mid grey tone
				glUniform3f(curS->getUniform("MatSpec"), 0.1f, 0.1f, 0.1f);    // specular: low highlights
				glUniform1f(curS->getUniform("MatShine"), 4.0f);               // low shininess for a matte finish
				break;
			case 10: // emissive star material (visible against dark backgrounds)
				glUniform3f(curS->getUniform("MatAmb"), 0.7f, 0.7f, 0.6f);   // Higher ambient for visibility in darkness
				glUniform3f(curS->getUniform("MatDif"), 1.0f, 1.0f, 0.8f);   // Bright yellow diffuse
				glUniform3f(curS->getUniform("MatSpec"), 0.5f, 0.5f, 0.4f);  // Higher specular
				glUniform1f(curS->getUniform("MatShine"), 15.0f);            // Moderate shininess
				break;
  		}
	}

	/* helper function to set model trasnforms */
  	void SetModel(vec3 trans, float rotY, float rotX, float sc, shared_ptr<Program> curS) {
  		mat4 Trans = glm::translate( glm::mat4(1.0f), trans);
  		mat4 RotX = glm::rotate( glm::mat4(1.0f), rotX, vec3(1, 0, 0));
  		mat4 RotY = glm::rotate( glm::mat4(1.0f), rotY, vec3(0, 1, 0));
  		mat4 ScaleS = glm::scale(glm::mat4(1.0f), vec3(sc));
  		mat4 ctm = Trans*RotX*RotY*ScaleS;
  		glUniformMatrix4fv(curS->getUniform("M"), 1, GL_FALSE, value_ptr(ctm));
  	}

	void setModel(std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack>M) {
		glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(M->topMatrix()));
   	}

	glm::mat4 computeShadowMatrix(const glm::mat4& modelMatrix, const glm::vec3& lightDir, float offsetFactor, float groundY) {
		glm::mat4 flatten = glm::scale(glm::mat4(1.0f), glm::vec3(1.0f, 0.0f, 1.0f));	// create flattening scalar matrix
		glm::mat4 flatModel = flatten * modelMatrix;									// flatten the model
		glm::vec3 projTranslation = glm::vec3(flatModel[3]);							// get translation of model
		projTranslation.y = groundY + 0.01f;											// set shadows y to just above ground height
		flatModel[3] = glm::vec4(projTranslation, 1.0f);								// set translation
		glm::vec3 lightHoriz = glm::normalize(glm::vec3(lightDir.x, 0.0f, lightDir.z)); // compute horizontal light direction
		glm::vec3 offset = -lightHoriz * offsetFactor;									// compute offset in oposite direction
		glm::mat4 offsetMat = glm::translate(glm::mat4(1.0f), offset);					// set translation with offset
		return offsetMat * flatModel;													// final matrix is a flat model with an offset
	}

	void render() {
		int width, height;
		string modelKey;

		glfwGetFramebufferSize(windowManager->getHandle(),
			&width, &height);										// get current frame buffer size
		glViewport(0, 0, width, height);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);			// clear framebuffer

		float aspect = width / (float)height;	

		// using the matrix stack from Lab 6
		// create the matrix stacks - please leave these alone for now
		auto Projection = make_shared<MatrixStack>();
		auto View = make_shared<MatrixStack>();
		auto Model = make_shared<MatrixStack>();

		Projection->pushMatrix();
		Projection->perspective(45.0f, aspect, 0.01f, 100.0f);
		View->pushMatrix();
		View->loadIdentity();
		
		// Compute delta time for the tour update
		float currentTime = glfwGetTime();
		float deltaTime = currentTime - lastTourUpdateTime;
		lastTourUpdateTime = currentTime;

		// compute the moon orbit position and light position
		float angle = currentTime * moonSpeed; // moon’s angular position
		float R = groundSize * multiMeshes["sphereWTex"].scale; // radius for moon’s orbit
		float moonRadius = sqrt(R * R - moonY * moonY);
		glm::vec3 moonCenter = (multiMeshes["spaceKit"].shapes[1]->max + multiMeshes["spaceKit"].shapes[1]->min) * 0.5f;
		glm::vec3 moonLightPos = glm::vec3(moonRadius * cos(angle), moonY, moonRadius * sin(angle));
		// FOR WHATEVER REASON THE BOUNDING BOX LIED AND I HAD TO DO THIS BY HAND
		moonCenter.z += 16.0f;
		moonCenter.x += 2.5f;
		glm::vec3 moonPos = moonLightPos - moonCenter;
		glm::vec3 lightDir = glm::normalize(moonLightPos); // direction of light for shadows

		// If tour mode is active, update the spline and set the eye position
		if (tourEnabled && cameraTour != nullptr) {
			cameraTour->update(deltaTime);
			// Set the camera's eye to the current position along the spline
			eye = cameraTour->getPosition();
			// Fix the lookAt point for the entire tour
			lookAt = tourLookAt;
			// If the spline is done, disable tour mode
			if (cameraTour->isDone()) {
				tourEnabled = false;
				delete cameraTour;
				cameraTour = nullptr;
			}
		} else {
			// Otherwise, update lookAt normally (if you are in manual mode)
			// Calculate lookAt if not initialized
			if (!mouseInitialized) {
				updateLookAt();
				mouseInitialized = true;
			} else {
				updateLookAt();
			}
		}

		// Use glm::lookAt to create view matrix
		glm::mat4 viewMatrix = glm::lookAt(eye, lookAt, glm::vec3(0.0f, 1.0f, 0.0f));
		View->multMatrix(viewMatrix);

		texProg->bind();											// switch to texture shader
		glUniformMatrix4fv(texProg->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
		glUniformMatrix4fv(texProg->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));
		glUniformMatrix4fv(texProg->getUniform("M"), 1, GL_FALSE, value_ptr(Model->topMatrix()));
		glUniform3f(texProg->getUniform("lightPos"), moonLightPos.x, moonLightPos.y, moonLightPos.z);
		glUniform1f(texProg->getUniform("MatShine"), 3.0);

		glDepthFunc(GL_LEQUAL);										// set the depth function to always draw the sphere!
		modelKey = "sphereWTex";									// sphere mesh for background
		texture1->bind(texProg->getUniform("Texture"));				// night sky texture
		glUniform1i(texProg->getUniform("flip"), 1);				// set flip normals flag
		if (multiMeshes.find(modelKey) != multiMeshes.end() && !multiMeshes[modelKey].shapes.empty()) {
			Model->pushMatrix();
			Model->loadIdentity();
			Model->scale(vec3(groundSize
				* multiMeshes[modelKey].scale));					// set scale factor
			setModel(texProg, Model);								// apply transforms
			multiMeshes[modelKey].shapes[0]->draw(texProg);			// draw the sphere
			Model->popMatrix();
		}
		glUniform1i(texProg->getUniform("flip"), 0);				// reset flip normals flag
		glDepthFunc(GL_LESS);										// set the depth test back to normal!

		drawGround(texProg);										// draw ground wrapper

		texProg->unbind();

		prog->bind(); // load the simple shader (no texture)
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
		glUniformMatrix4fv(prog->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));
		//glUniform3f(prog->getUniform("lightPos"), -2.0 + lightTrans, 2.0, 2.0);
		// update light position uniform in shaders so that the moon is the light source
		glUniform3f(prog->getUniform("lightPos"), moonLightPos.x, moonLightPos.y, moonLightPos.z);
		glUniform3f(prog->getUniform("lightIntensity"), 0.6, 0.6, 0.65);
		// set a low level ambient light
		glUniform3f(prog->getUniform("ambientLight"), 0.4f, 0.4f, 0.6f);
		glUniform1f(prog->getUniform("uTime"), glfwGetTime());
		glUniform1f(prog->getUniform("windStrength"), 0.02f);
		glUniform1f(prog->getUniform("windFrequency"), 2.0f);
		glUniform1i(prog->getUniform("applyWind"), 0); // false by default
		glUniform1i(prog->getUniform("useFixedLight"), 0); // false by default
		glUniform3f(prog->getUniform("fixedLightPos"), 0, moonY, 0);
		glUniform1i(prog->getUniform("ignoreLight"), 0); // defualt ignore lighting is off


		glUniform1i(prog->getUniform("applyWind"), 1);				// set wave animation flag to true
		modelKey = "foliageCollection";								// small plant mesh collection
		if (foliagePlacements.find(modelKey) != foliagePlacements.end() &&
			multiMeshes.find(modelKey) != multiMeshes.end()) {
			MultiMesh& foliage = multiMeshes[modelKey];				// get meshes
			vector<vector<glm::mat4>>& placements =
				foliagePlacements[modelKey];						// get precomputed placements
			for (size_t i = 0; i < foliage.shapes.size(); i++) {	// for each plant type (sub-mesh)
				for (const auto& modelMatrix : placements[i]) {		// for each instance of this plant type
					SetMaterial(prog, 2);							// apply plant material
					glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(modelMatrix));
					foliage.shapes[i]->draw(prog);					// draw plant mesh
				}
			}
		}
		glUniform1i(prog->getUniform("applyWind"), 0);				// set wave animation flag to false

		// set the flag so the moon uses a global light
		glUniform1i(prog->getUniform("useFixedLight"), 1);

		glDepthFunc(GL_LEQUAL);  // allow equal depth values to pass the test

		modelKey = "spaceKit";											// moon mesh
		if (multiMeshes.find(modelKey) != multiMeshes.end() && !multiMeshes[modelKey].shapes.empty()) {
			Model->pushMatrix();
			Model->loadIdentity();
			Model->translate(moonPos);							// translate to orbit position
			Model->scale(vec3(2.0f));								// set scale factor
			setModel(prog, Model);									// apply transforms
			SetMaterial(prog, 1);									// apply moon material
			multiMeshes[modelKey].shapes[0]->draw(prog);			// draw moon
			multiMeshes[modelKey].shapes[1]->draw(prog);			// draw moon
			Model->popMatrix();
		}

		glDepthFunc(GL_LESS); // restore normal depth testing

		updateStarAnimation(currentTime); // update star animations before they get drawn
		modelKey = "star"; // star mesh

		// disable depth writing but keep depth testing
		glDepthMask(GL_FALSE);

		for (int i = 0; i < starPlacements.size(); i++) {
			// if the star has faded completely dont draw it
			if (starAlphas[i] < 0.01f) continue;

			glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(starPlacements[i])); // set model matrix
			SetMaterial(prog, 3); // set star material

			// set alpha for this star
			glUniform4f(prog->getUniform("starColor"), 1.0f, 1.0f, 0.8f, starAlphas[i]);

			// enable blending with additive mode
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);

			multiMeshes[modelKey].shapes[0]->draw(prog); // draw star

			// disable blending
			glDisable(GL_BLEND);
		}

		// restore OpenGL state
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);

		// reset the flag signaling global light for other objects
		glUniform1i(prog->getUniform("useFixedLight"), 0);

		modelKey = "natureCollection";  // nature mesh collection for trees
		if (treePlacements.find(modelKey) != treePlacements.end() &&
			multiMeshes.find(modelKey) != multiMeshes.end()) {
			MultiMesh& nature = multiMeshes[modelKey];
			vector<vector<glm::mat4>>& placements = treePlacements[modelKey];
			for (size_t i = 0; i < placements.size(); i++) {
				for (const auto& modelMatrix : placements[i]) {
					glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(modelMatrix));
					if (treeIndices[i] == 11) {						// indice of tree
						SetMaterial(prog, 5);						// set material to brown tree color
						nature.shapes[treeIndices[i]]->draw(prog);	// draw trunk

						glUniform1i(prog->getUniform("ignoreLight"), 1); // set to ignore lighting
						glm::mat4 shadowMat = computeShadowMatrix(modelMatrix, lightDir, 0.2f, -1.25f);
						glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(shadowMat));
						SetMaterial(prog, 8);						// shadow material
						nature.shapes[treeIndices[i]]->draw(prog);	// draw trunk shadow
						glUniform1i(prog->getUniform("ignoreLight"), 0); // reset ignore lighting

						glUniform1i(prog->getUniform("applyWind"), 1);
						glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(modelMatrix));
						SetMaterial(prog, 6);						// set material to leaves color
						nature.shapes[treeIndices[i]
							+ 1]->draw(prog);						// draw leaves

						glUniform1i(prog->getUniform("ignoreLight"), 1); // set to ignore lighting
						shadowMat = computeShadowMatrix(modelMatrix, lightDir, 0.2f, -1.25f);
						glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(shadowMat));
						SetMaterial(prog, 8);						// shadow material
						nature.shapes[treeIndices[i] + 1]->draw(prog);	// draw leaves shadow
						glUniform1i(prog->getUniform("ignoreLight"), 0); // reset ignore lighting

						glUniform1i(prog->getUniform("applyWind"), 0);
					} else {										// logs and stumps drawn normally
						SetMaterial(prog, 5);						// set material to brown tree color
						nature.shapes[treeIndices[i]]->draw(prog);	// draw logs and stumps
					}
				}
			}
		}

		for (auto& entry : rockPlacements) {
			string rockKey = entry.first;
			if (multiMeshes.find(rockKey) == multiMeshes.end()) continue;
			for (const auto& modelMatrix : entry.second) {
				glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, glm::value_ptr(modelMatrix));
				SetMaterial(prog, 0);								// stone like material
				multiMeshes[rockKey].shapes[0]->draw(prog);			// draw rock

				glUniform1i(prog->getUniform("ignoreLight"), 1); // set to ignore lighting
				glm::mat4 shadowMat = computeShadowMatrix(modelMatrix, lightDir, 0.15f, -1.25f);
				glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(shadowMat));
				SetMaterial(prog, 9);								// shadow material
				multiMeshes[rockKey].shapes[0]->draw(prog);			// draw shadows for rocks
				glUniform1i(prog->getUniform("ignoreLight"), 0); // reset ignore lighting
			}
		}

		modelKey = "moonDoor";									// door mesh
		if (multiMeshes.find(modelKey) != multiMeshes.end() && !multiMeshes[modelKey].shapes.empty()) {
			Model->pushMatrix();
				Model->loadIdentity();
				Model->translate(doorPos);								// translate to the door position
				Model->scale(2.0f * glm::vec3(multiMeshes[modelKey].scale));	// scale the door to scale factor
				setModel(prog, Model);									// Set the model matrix
				SetMaterial(prog, 7);									// Set door mat
				multiMeshes[modelKey].shapes[0]->draw(prog);			// draw door

				Model->pushMatrix();
					glUniform1i(prog->getUniform("ignoreLight"), 1); // set to ignore lighting
					glm::mat4 shadowMat = computeShadowMatrix(Model->topMatrix(), lightDir, 0.15f, -1.25f);
					glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(shadowMat));
					SetMaterial(prog, 9);								// shadow material
					multiMeshes[modelKey].shapes[0]->draw(prog);			// draw shadows for rocks
					glUniform1i(prog->getUniform("ignoreLight"), 0); // reset ignore lighting
				Model->popMatrix();
			Model->popMatrix();
		}
		prog->unbind();


		// Pass the emitter position to the particle system
		doorParticleSystem->setEmitter(doorPos);

		// Draw
		partProg->bind();
		texture2->bind(partProg->getUniform("alphaTexture"));
		CHECKED_GL_CALL(glUniformMatrix4fv(partProg->getUniform("P"), 1, GL_FALSE, value_ptr(P->topMatrix())));
		CHECKED_GL_CALL(glUniformMatrix4fv(partProg->getUniform("V"), 1, GL_FALSE, value_ptr(V->topMatrix())));
		CHECKED_GL_CALL(glUniformMatrix4fv(partProg->getUniform("M"), 1, GL_FALSE, value_ptr(M->topMatrix())));

		CHECKED_GL_CALL(glUniform3f(partProg->getUniform("pColor"), 0.9, 0.7, 0.7));

		doorParticleSystem->drawMe(partProg);
		doorParticleSystem->update();

		partProg->unbind();

		// pop matrix stacks
		Projection->popMatrix();
		View->popMatrix();
	}
};

void mouseMoveCallbackWrapper(GLFWwindow* window, double xpos, double ypos) {
	Application* app = (Application*)glfwGetWindowUserPointer(window);
	app->mouseMoveCallback(window, xpos, ypos);
}

void scrollCallbackWrapper(GLFWwindow* window, double deltaX, double deltaY) {
	Application* app = (Application*)glfwGetWindowUserPointer(window);
	app->scrollCallback(window, deltaX, deltaY);
}

void mouseButtonCallbackWrapper(GLFWwindow* window, int button, int action, int mods) {
	Application* app = (Application*)glfwGetWindowUserPointer(window);
	app->mouseCallback(window, button, action, mods);
}

int main(int argc, char *argv[]) {
	std::string resourceDir = "../resources";	// Where the resources are loaded from

	if (argc >= 2) {
		resourceDir = argv[1];
	}

	Application *application = new Application();
	WindowManager *windowManager = new WindowManager();
	windowManager->init(1920, 1080);
	windowManager->setEventCallbacks(application);
	application->windowManager = windowManager;

	glfwSetWindowUserPointer(windowManager->getHandle(), application);
	glfwSetCursorPosCallback(windowManager->getHandle(), mouseMoveCallbackWrapper);
	glfwSetScrollCallback(windowManager->getHandle(), scrollCallbackWrapper);
	glfwSetMouseButtonCallback(windowManager->getHandle(), mouseButtonCallbackWrapper);

	application->init(resourceDir);
	application->initGeom(resourceDir);
	application->initTex(resourceDir);
	application->initGround();										// code to load in the ground plane (CPU defined data passed to GPU)
	application->initFoliage();
	application->initStars();
	application->initNature();
	application->initRocks();

	while (! glfwWindowShouldClose(windowManager->getHandle())) {	// Loop until the user closes the window.
		application->render();										// Render scene
		glfwSwapBuffers(windowManager->getHandle());				// Swap front and back buffers
		glfwPollEvents();											// Poll for and process events
	}
	windowManager->shutdown();										// Quit program
	return 0;
}