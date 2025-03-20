/* Final Project
 * Max Eyrich
 * Built off P4 base code - Cal Poly Z. Wood + S. Sueda + I. Dunn
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
#include "FBXModel.h"
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;
using namespace std;
using namespace glm;

class Application : public EventCallbacks {
public:

	WindowManager * windowManager = nullptr;

	// shader programs
	std::shared_ptr<Program> prog, texProg, waterProg;

	// geometry
	std::vector<std::shared_ptr<Shape>> meshes;

	struct MultiMesh {
		vector<shared_ptr<Shape>> shapes; // Each sub-mesh in the object
		float scale; // Uniform scale for multimesh objects
	};

	unordered_map<string, MultiMesh> multiMeshes; // Store multi-mesh objects

	// wolf data
	std::shared_ptr<Program> animShader;
	FBXModel* wolfModel;
	float lastTime;
	float wolfSpeed;

	// containers for static placements
	std::vector<glm::mat4> starPlacements;
	unordered_map<string, vector<vector<glm::mat4>>> foliagePlacements, treePlacements;
	unordered_map<string, vector<glm::mat4>> rockPlacements;

	// global data for ground plane - direct load constant defined CPU data to GPU (not obj)
	GLuint GrndBuffObj, GrndNorBuffObj, GrndTexBuffObj, GIndxBuffObj;
	int g_GiboLen;
	// ground VAO
	GLuint GroundVertexArrayID;

	// global data for water
	GLuint WaterVAO, WaterVBO, WaterEBO, WaterIBO;
	int waterIndexCount;

	// textures
	shared_ptr<Texture> groundTex; // Ground texture
	shared_ptr<Texture> skyTex; // Sky texture
	shared_ptr<Texture> doorTex; // brick texture
	shared_ptr<Texture> basinTex; // Sandy texture for the basin

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
	float tourDuration = 15.0f;            // Duration (in seconds) for the tour
	float lastTourUpdateTime = 0.0f;       // For computing delta time in tour mode

	float groundSize = 20.0f;
	float sphereRadius = groundSize / 2.0f;
	float sphereEdge;

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
	float lightTrans = 0;
	float moonSpeed = 0.15f;						// angular speed (rad/s)
	float moonY = sphereRadius * 0.66f;				// orbit height

	// door lighting
	float doorLightIntensity = 1.0f;

	// particle system
	particleSys* doorParticleSystem = nullptr;
	std::shared_ptr<Program> partProg; // Shader program for particles
	std::shared_ptr<Texture> particleTexture; // Texture for particles

	void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
	{
		if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) { // exit the program when 'ESC' is pressed
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
		if (key == GLFW_KEY_G && action == GLFW_PRESS) { // toggle cinematic tour when 'g' is pressed
			if (!tourEnabled) {	// start tour
				tourLookAt = glm::vec3(0.0f, 0.0f, 0.0f); // fix the look-at point
				glm::vec3 start = vec3(0, 7.0f, sphereEdge); // start tour above doorPos
				glm::vec3 control = start + glm::vec3(sphereEdge, 3.0f, 0.0f); // choose a control point to create a nice arc
				glm::vec3 end = glm::vec3(0.0f, 7.0f, -sphereEdge); // end point (destination) for camera to travel to
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
		if (key == GLFW_KEY_N && action == GLFW_PRESS) { // Switch to next wolf animation
			int nextAnim = (wolfModel->getCurrentAnimation() + 1) % wolfModel->getAnimationCount();
			wolfModel->setAnimation(nextAnim);
		}
		if (key == GLFW_KEY_P && action == GLFW_PRESS) { // Switch to previous animation
			int prevAnim = wolfModel->getCurrentAnimation() - 1;
			if (prevAnim < 0) prevAnim = wolfModel->getAnimationCount() - 1;
			wolfModel->setAnimation(prevAnim);
		}
		if (key == GLFW_KEY_0 && action == GLFW_PRESS && mods == GLFW_MOD_CONTROL) { // Print out all animation names
			std::cout << "Available animations:" << std::endl;
			for (int i = 0; i < wolfModel->getAnimationCount(); i++) {
				std::cout << "  " << i << ": " << wolfModel->getAnimationName(i) << std::endl;
			}
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

		texProg = make_shared<Program>();
		texProg->setVerbose(true);
		texProg->setShaderNames(resourceDirectory + "/tex_vert.glsl", resourceDirectory + "/tex_frag.glsl");
		texProg->init();
		texProg->addUniform("P");
		texProg->addUniform("V");
		texProg->addUniform("M");
		texProg->addUniform("flip");
		texProg->addUniform("ground");
		texProg->addUniform("basin");
		texProg->addUniform("Texture");
		texProg->addUniform("MatShine");
		texProg->addUniform("lightPos");
		texProg->addAttribute("vertPos");
		texProg->addAttribute("vertNor");
		texProg->addAttribute("vertTex");

		waterProg = make_shared<Program>();
		waterProg->setVerbose(true);
		waterProg->setShaderNames(resourceDirectory + "/water_vert.glsl", resourceDirectory + "/water_frag.glsl");
		waterProg->init();
		waterProg->addUniform("P");
		waterProg->addUniform("V");
		waterProg->addUniform("M");
		waterProg->addUniform("uTime");
		waterProg->addUniform("lightPos");
		waterProg->addUniform("viewPos");
		waterProg->addAttribute("aPos");

		// Initialize animation shader
		animShader = make_shared<Program>();
		animShader->setVerbose(true);
		animShader->setShaderNames(resourceDirectory + "/anim_vert.glsl", resourceDirectory + "/anim_frag.glsl");
		animShader->init();
		animShader->addUniform("P");
		animShader->addUniform("V");
		animShader->addUniform("M");
		animShader->addUniform("MatAmb");
		animShader->addUniform("MatDif");
		animShader->addUniform("MatSpec");
		animShader->addUniform("MatShine");
		animShader->addUniform("lightPos");

		// Add uniforms for bone transforms
		for (int i = 0; i < 100; i++) {
			animShader->addUniform("boneTransforms[" + std::to_string(i) + "]");
		}

		animShader->addAttribute("vertPos");
		animShader->addAttribute("vertNor");
		animShader->addAttribute("vertTex");
		animShader->addAttribute("boneIDs");
		animShader->addAttribute("weights");
		animShader->addUniform("hasTexture");
		animShader->addUniform("textureSampler");

		// Create the wolf model
		wolfModel = new FBXModel();
		lastTime = glfwGetTime();
	}

	void initGeom(const std::string& resourceDirectory) {
		// mesh objects to target
		vector<string> multiMeshFiles = {
			"foliageCollection.obj",
			"moon.obj",
			"natureCollection.obj",
			"rock1.obj",
			"rock2.obj",
			"rock3.obj",
			"rock4.obj",
			"rock5.obj",
			"sphereWTex.obj",
			"star.obj",
			"spaceKit.obj",
			"moonDoor.obj",
			"cube.obj"
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

	void initParticleSystem(const std::string& resourceDirectory) {
		// Initialize the particle shader program
		partProg = make_shared<Program>();
		partProg->setVerbose(true);
		partProg->setShaderNames(
			resourceDirectory + "/lab10_vert.glsl",
			resourceDirectory + "/lab10_frag.glsl");
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

		// Load the particle texture
		particleTexture = make_shared<Texture>();
		particleTexture->setFilename(resourceDirectory + "/alpha.bmp");
		particleTexture->init();
		particleTexture->setUnit(3); // Use texture unit 3
		particleTexture->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);

		// Create the particle system (initially at origin)
		doorParticleSystem = new particleSys(vec3(0, 0, 0));
		doorParticleSystem->gpuSetup();
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
			float R = sphereRadius * multiMeshes[sphereKey].scale * 0.99f; // small inner radius offset
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
		float rInner = 5.5f;
		float rOuter = sphereRadius * multiMeshes["sphereWTex"].scale * 2.0f;

		// setup random number generators
		std::default_random_engine rng(std::random_device{}());
		std::uniform_int_distribution<int> countDist(100, 300);						// random range of plants per type
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
		vector<vector<glm::mat4>> placements;   // prepare a container for the placements
		placements.resize(3);                   // one vector per tree type

		// define inner and outer radii for annular region
		float rInner = 4.0f * sphereRadius * multiMeshes["sphereWTex"].scale / 3.0f;
		float rOuter = 2.0f * sphereRadius * multiMeshes["sphereWTex"].scale - 5.0f;

		// Set up random generators
		std::default_random_engine rng(std::random_device{}());
		std::uniform_real_distribution<float> rDist(rInner, rOuter);
		std::uniform_real_distribution<float> thetaDist(0.0f, glm::two_pi<float>());
		std::uniform_int_distribution<int> countDistT(50, 70);                     // random number of trees
		std::uniform_int_distribution<int> countDistLS(10, 20);                     // random number of stumps and logs
		std::uniform_real_distribution<float> rotDist(0.0f, glm::two_pi<float>());  // random rotation

		// Define door position and the path that should be kept clear
		glm::vec3 doorPos(0, -1.25f, sphereEdge); // Door position
		glm::vec3 centerPos(0, -1.25f, 0); // Center of the scene

		// Calculate door-to-center direction vector
		glm::vec3 doorToCenter = glm::vec3(0, 0, -1.0f);

		float corridorWidth = 0.0f; // Width of clear path

		for (int i = 0; i < 3; i++) { // range of placements per type
			int instanceCount = i < 1 ? countDistT(rng) : countDistLS(rng);
			// Vector to store valid placements for this tree type
			std::vector<glm::mat4> validPlacements;

			// Try to place more trees than needed, then filter out those in the corridor
			int attemptsLimit = instanceCount * 2; // Generate twice as many candidates

			for (int j = 0; j < attemptsLimit && validPlacements.size() < instanceCount; j++) {
				// Choose a random radius and angle within the annular region
				float r = rDist(rng);
				float theta = thetaDist(rng);
				float x = r * cos(theta);
				float z = r * sin(theta);
				float y = -1.25f; // Fix the vertical position

				// Create potential position
				glm::vec3 treePos(x, y, z);

				// Check if the tree would block the path from door to center
				// Calculate distance from point to line (door-to-center)
				glm::vec3 doorToTree = treePos - doorPos;
				float projectionLength = glm::dot(doorToTree, doorToCenter);

				// Only consider trees that are between the door and center (not behind the door)
				bool inPathDirection = projectionLength > 0 && projectionLength < sphereEdge;

				// Calculate perpendicular distance from tree to the door-center line
				glm::vec3 projection = doorPos + doorToCenter * projectionLength;
				float distanceToPath = glm::length(treePos - projection);

				// Skip this position if the tree would be in the corridor
				if (inPathDirection && distanceToPath < corridorWidth) {
					continue;  // Skip this position
				}

				// Position is valid, create the model matrix
				float angle = rotDist(rng);                                        // Random rotation about Y
				glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z));  // Build model matrix and set translation
				model = glm::rotate(model, angle, glm::vec3(0, 1, 0));             // Set random rotation
				model = glm::scale(model, 5.0f * glm::vec3(1, 1, 1));              // Set scale

				validPlacements.push_back(model);
			}

			// Store the valid placements (might be less than requested if many were filtered out)
			placements[i] = validPlacements;
		}

		treePlacements[natureKey] = placements;
	}

	void initRocks() {
		vector<string> rockKeys = { "rock1", "rock2", "rock3", "rock4", "rock5" };

		// define inner and outer radii for annular region
		float rInner = 5.25f;
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

	void initRingRocks() {
		vector<string> rockKeys = { "rock1", "rock2", "rock3", "rock4", "rock5" };

		// Define parameters for the rock ring
		float basinRadius = 5.0f; // Same as holeRadius in initWaterSurface
		float ringRadius = basinRadius - 0.1f; // Slightly smaller than the water basin
		int numRocks = 100; // Number of rocks to place around the ring

		// Set up random generators for variation
		std::default_random_engine rng(std::random_device{}());
		std::uniform_real_distribution<float> radiusVariation(-0.2f, 0.2f); // Small variation in radius
		std::uniform_real_distribution<float> rotDist(0.0f, glm::two_pi<float>()); // Random rotation
		std::uniform_real_distribution<float> scaleDist(0.8f, 1.3f); // Random scaling
		std::uniform_real_distribution<float> heightDist(-0.1f, 0.1f); // Small height variation
		std::uniform_int_distribution<int> rockTypeDist(0, rockKeys.size() - 1); // Random rock type

		for (int i = 0; i < numRocks; i++) {
			// Choose a random rock type from the available ones
			string rockKey = rockKeys[rockTypeDist(rng)];
			if (multiMeshes.find(rockKey) == multiMeshes.end())
				continue;

			auto rockShape = multiMeshes[rockKey].shapes[0];
			float rockYOffset = -1.5f - rockShape->min.y; // Y offset

			// Calculate position around the circle
			float angle = (float)i / numRocks * glm::two_pi<float>();
			float r = ringRadius + radiusVariation(rng); // Add small variation to radius
			float x = r * cos(angle);
			float z = r * sin(angle);
			float y = rockYOffset + heightDist(rng); // Add small height variation

			// Create model matrix with translation, rotation and scaling
			glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(x, y, z));
			model = glm::rotate(model, rotDist(rng), glm::vec3(0, 1, 0)); // Random Y rotation

			// Add small random rotations on X and Z for more natural look
			model = glm::rotate(model, radiusVariation(rng) * 0.3f, glm::vec3(1, 0, 0));
			model = glm::rotate(model, radiusVariation(rng) * 0.3f, glm::vec3(0, 0, 1));

			// Apply scaling - smaller scale for a decorative ring
			float scale = scaleDist(rng) * 0.7f; // Smaller scale than the other scattered rocks
			model = glm::scale(model, glm::vec3(scale));

			// Add this placement to rock placements
			if (rockPlacements.find(rockKey) == rockPlacements.end()) {
				rockPlacements[rockKey] = vector<glm::mat4>();
			}
			rockPlacements[rockKey].push_back(model);
		}
		std::cout << "Added ring of rocks around water basin" << std::endl;
	}

	// code to load in the textures
	void initTex(const std::string& resourceDirectory) {
		groundTex = make_shared<Texture>();
		groundTex->setFilename(resourceDirectory + "/moss_ground_2k.jpg");
		groundTex->init();
		groundTex->setUnit(0);
		groundTex->setWrapModes(GL_REPEAT, GL_REPEAT);

		skyTex = make_shared<Texture>();
		skyTex->setFilename(resourceDirectory + "/starSky.jpg");
		skyTex->init();
		skyTex->setUnit(1);
		skyTex->setWrapModes(GL_REPEAT, GL_REPEAT);

		doorTex = make_shared<Texture>();
		doorTex->setFilename(resourceDirectory + "/stone_bricks_2k.jpg");
		doorTex->init();
		doorTex->setUnit(2);
		doorTex->setWrapModes(GL_REPEAT, GL_REPEAT);

		basinTex = make_shared<Texture>();
		basinTex->setFilename(resourceDirectory + "/sand_2k.jpg");
		basinTex->init();
		basinTex->setUnit(4);
		basinTex->setWrapModes(GL_REPEAT, GL_REPEAT);
	}

	void initWaterSurface(float radius) {
		// Create a flat circular plane for the water
		std::vector<float> positions;
		std::vector<unsigned int> indices;

		int segments = 32;

		// Center vertex
		positions.push_back(0.0f); // x
		positions.push_back(-1.35f); // y (slightly below ground level)
		positions.push_back(0.0f); // z

		// Edge vertices
		for (int i = 0; i <= segments; i++) {
			float angle = 2.0f * glm::pi<float>() * i / segments;
			float x = radius * cos(angle);
			float z = radius * sin(angle);

			positions.push_back(x);
			positions.push_back(-1.25f); // Same level as hole in ground
			positions.push_back(z);
		}

		// Create triangles with center as one vertex
		for (int i = 1; i <= segments; i++) {
			indices.push_back(0); // Center
			indices.push_back(i);
			indices.push_back(i + 1 > segments ? 1 : i + 1);
		}

		// Create VAO and VBO
		glGenVertexArrays(1, &WaterVAO);
		glBindVertexArray(WaterVAO);

		glGenBuffers(1, &WaterVBO);
		glBindBuffer(GL_ARRAY_BUFFER, WaterVBO);
		glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), positions.data(), GL_STATIC_DRAW);

		glGenBuffers(1, &WaterIBO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, WaterIBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

		waterIndexCount = indices.size();

		// Position attribute
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

		glBindVertexArray(0);
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

    // code to draw the ground plane
    void drawGround(shared_ptr<Program> curS) {
		curS->bind();
		glBindVertexArray(GroundVertexArrayID);
		glUniform1i(curS->getUniform("ground"), 5); // enable the hole
		glUniform1i(curS->getUniform("basin"), 0);  // normal drawing mode
		glUniform1f(texProg->getUniform("MatShine"), 20.0);
		groundTex->bind(curS->getUniform("Texture"));
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
		glUniform1i(curS->getUniform("ground"), 0); // disable the hole

		// Draw the basin (half-sphere) below the hole
		string sphereKey = "sphereWTex";
		if (multiMeshes.find(sphereKey) != multiMeshes.end() && !multiMeshes[sphereKey].shapes.empty()) {
			// Switch to the sandy texture
			basinTex->bind(curS->getUniform("Texture"));

			// Position and scale the basin
			float basinDepth = 2.0f;
			float holeRadius = 5.0f;

			mat4 basinModel = glm::translate(glm::mat4(1.0f), vec3(0, -1.25f, 0));

			// Scale the sphere to match the hole radius
			float scaleXZ = holeRadius / multiMeshes[sphereKey].shapes[0]->max.x;
			basinModel = glm::scale(basinModel, vec3(scaleXZ, scaleXZ, scaleXZ));

			// Set model matrix for the basin
			glUniformMatrix4fv(curS->getUniform("M"), 1, GL_FALSE, value_ptr(basinModel));

			glUniform1i(curS->getUniform("basin"), 1);	// Enable basin drawing mode (clip above ground)
			glUniform1i(curS->getUniform("ground"), 0); // No hole for basin
			glUniform1i(curS->getUniform("flip"), 1);   // Flip normals

			multiMeshes[sphereKey].shapes[0]->draw(curS); // Draw the basin

			// Reset uniforms
			glUniform1i(curS->getUniform("basin"), 0);
			glUniform1i(curS->getUniform("flip"), 0);
		}

  		curS->unbind(); // unbind shader
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
		glm::vec3 moonPos = moonLightPos - moonCenter;

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
		glUniform1f(texProg->getUniform("MatShine"), 1.0);
		glUniform1i(texProg->getUniform("ground"), 0); // disable the hole
		glUniform1i(texProg->getUniform("basin"), 0);

		glDepthFunc(GL_LEQUAL);										// set the depth function to always draw the sphere!
		modelKey = "sphereWTex";									// sphere mesh for background
		skyTex->bind(texProg->getUniform("Texture"));				// night sky texture
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

		prog->bind();												// load the simple shader (no texture)
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

		glm::vec3 lightDir = glm::normalize(moonLightPos);			// direction of light for shadows

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

		prog->unbind();

		// Render the moonDoor
		modelKey = "moonDoor";  // moonDoor mesh key
		if (multiMeshes.find(modelKey) != multiMeshes.end() && !multiMeshes[modelKey].shapes.empty()) {
			//===========
			// draw door
			//===========
			Model->pushMatrix();
				Model->loadIdentity();

				// Position the door on the edge of the scene
				glm::vec3 doorPos(0, -1.25f, sphereEdge);
				Model->translate(doorPos);

				Model->scale(vec3(multiMeshes[modelKey].scale * 2.0f)); // Scale appropriate to the scene
				
				texProg->bind(); // Use texture shader since this will be textured

				// Set up uniforms for the texture program
				glUniformMatrix4fv(texProg->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
				glUniformMatrix4fv(texProg->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));
				glUniform3f(texProg->getUniform("lightPos"), moonLightPos.x, moonLightPos.y, moonLightPos.z);
				glUniform1f(texProg->getUniform("MatShine"), 10.0);
				glUniform1i(texProg->getUniform("ground"), 0); // disable the hole
				glUniform1i(texProg->getUniform("basin"), 0);

				doorTex->bind(texProg->getUniform("Texture")); // Create texture for the door

				setModel(texProg, Model); // Apply the model matrix

				multiMeshes[modelKey].shapes[0]->draw(texProg); // Draw the moonDoor

				texProg->unbind(); // unbind texture shader

				//==============================
				// draw the illuminated doorway
				//==============================
				Model->pushMatrix();
					prog->bind(); // Use the non-textured shader
					float doorY = (multiMeshes["moonDoor"].shapes[0]->max.y + multiMeshes["moonDoor"].shapes[0]->min.y) * 0.5f;
					Model->translate(vec3(0.0f, doorY, 0.0f)); // Slightly in front of the door
					Model->scale(vec3(2.0f, 4.0f, 0.01f)); // Scale the light to fit the doorway opening

					// Apply animation with limited opacity range (30-50%)
					float minOpacity = 0.3f; // Minimum opacity (30%)
					float maxOpacity = 0.5f; // Maximum opacity (50%)
					float opacityRange = maxOpacity - minOpacity;
					float baseOpacity = minOpacity + (opacityRange / 2.0f); // Center point of the range
					float fluctuation = opacityRange / 2.0f; // Half of the range for fluctuation
					doorLightIntensity = baseOpacity + fluctuation * sin(currentTime * 0.3f);

					glUniformMatrix4fv(prog->getUniform("M"), 1, GL_FALSE, value_ptr(Model->topMatrix()));

					// Set up color with constant RGB but varying alpha
					glUniform3f(prog->getUniform("MatAmb"), 0.6f, 0.2f, 0.0f);     // deep orange-red ambient
					glUniform3f(prog->getUniform("MatDif"), 1.0f, 0.5f, 0.0f);     // bright orange diffuse
					glUniform3f(prog->getUniform("MatSpec"), 1.0f, 0.7f, 0.3f);    // golden yellow specular
					glUniform1f(prog->getUniform("MatShine"), 0.0f);               // no shininess for a pure glow

					// Enable blending for the glow effect with alpha controlling opacity
					glEnable(GL_BLEND);
					glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Standard alpha blending

					// Set the alpha value in the starColor uniform (assuming this uniform takes an RGBA)
					glUniform4f(prog->getUniform("starColor"), 1.0f, 0.6f, 0.0f, doorLightIntensity);

					multiMeshes["cube"].shapes[0]->draw(prog); // Draw the light plane (cube)

					glDisable(GL_BLEND); // Restore OpenGL state
					prog->unbind(); // unbind normal shader

				Model->popMatrix(); // Pop the light plane transforms

				//==============================
				// Draw the particles
				//==============================
				Model->pushMatrix();
					// Move particles to the center of the doorway
					glm::vec3 doorCenter = (multiMeshes["moonDoor"].shapes[0]->min + multiMeshes["moonDoor"].shapes[0]->max) * 0.5f;
					Model->translate(doorCenter);
					// Maybe slightly forward from the door plane
					Model->translate(vec3(0.0f, 0.0f, 0.1f));

					// Set the emitter position
					doorParticleSystem->setEmitter(vec3(0, 0, 0)); // Local space origin - transforms handled by matrix

					// Set the camera matrix
					doorParticleSystem->setCamera(View->topMatrix());

					// Draw the particles
					partProg->bind();
					particleTexture->bind(partProg->getUniform("alphaTexture"));
					CHECKED_GL_CALL(glUniformMatrix4fv(partProg->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix())));
					CHECKED_GL_CALL(glUniformMatrix4fv(partProg->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix())));
					CHECKED_GL_CALL(glUniformMatrix4fv(partProg->getUniform("M"), 1, GL_FALSE, value_ptr(Model->topMatrix())));

					// Set a warm color for the particles to match the door's glow
					CHECKED_GL_CALL(glUniform3f(partProg->getUniform("pColor"), 1.0f, 0.7f, 0.2f));

					doorParticleSystem->drawMe(partProg);
					doorParticleSystem->update();

					partProg->unbind();
				Model->popMatrix(); // Pop the particle transforms
			Model->popMatrix(); // pop the door transforms
		}

		//============
		// draw water
		//============
		waterProg->bind();
		glUniformMatrix4fv(waterProg->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
		glUniformMatrix4fv(waterProg->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));

		Model->pushMatrix();
			Model->loadIdentity();
			// Position at the center of the scene, slightly above ground
			Model->translate(vec3(0.0f, -0.24f, 0.0f));
			glUniformMatrix4fv(waterProg->getUniform("M"), 1, GL_FALSE, value_ptr(Model->topMatrix()));

			// Pass time for animation
			glUniform1f(waterProg->getUniform("uTime"), currentTime);
			// Pass light position and view position for reflections
			glUniform3f(waterProg->getUniform("lightPos"), moonLightPos.x, moonLightPos.y, moonLightPos.z);
			glUniform3f(waterProg->getUniform("viewPos"), eye.x, eye.y, eye.z);

			// Enable blending for transparency
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			// Draw the water surface
			glBindVertexArray(WaterVAO);
			glDrawElements(GL_TRIANGLES, waterIndexCount, GL_UNSIGNED_INT, 0);
			glBindVertexArray(0);

			// Disable blending
			glDisable(GL_BLEND);

		Model->popMatrix();
		waterProg->unbind();

		//===========
		// draw wolf
		//===========
		// Update wolf animation
		float wolfDeltaTime = currentTime - lastTime;
		lastTime = currentTime;

		wolfModel->update(wolfDeltaTime);

		// Update wolf position (walking in a circle)
		float radius = 6.0f;
		float wolfAngle = currentTime * wolfSpeed;

		// Calculate current position on the circle
		glm::vec3 wolfPosition;
		wolfPosition.x = radius * cos(wolfAngle);
		wolfPosition.y = -1.25f; // Ground level
		wolfPosition.z = radius * sin(wolfAngle);

		// Calculate direction of movement (tangent to the circle)
		glm::vec3 direction(-sin(wolfAngle), 0.0f, cos(wolfAngle));
		direction = glm::normalize(direction);

		// Set wolf position
		wolfModel->setPosition(wolfPosition);

		// Calculate the rotation angle to face the direction of movement
		// assumes the wolf's forward direction is +Z
		float rotationAngle = atan2(direction.x, direction.z) + glm::pi<float>();

		// Set wolf rotation
		wolfModel->setRotation(rotationAngle);

		// Draw the wolf
		animShader->bind();
		glUniformMatrix4fv(animShader->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
		glUniformMatrix4fv(animShader->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));
		glUniform3f(animShader->getUniform("MatAmb"), 0.1f, 0.1f, 0.1f);
		glUniform3f(animShader->getUniform("MatDif"), 0.8f, 0.8f, 0.8f);
		glUniform3f(animShader->getUniform("MatSpec"), 0.3f, 0.3f, 0.3f);
		glUniform1f(animShader->getUniform("MatShine"), 32.0f);
		glUniform3f(animShader->getUniform("lightPos"), moonLightPos.x, moonLightPos.y, moonLightPos.z);

		// Important: Initialize the texture uniforms
		glUniform1i(animShader->getUniform("hasTexture"), 0); // Will be set to 1 inside draw method if textures exist
		glUniform1i(animShader->getUniform("textureSampler"), 0); // Use texture unit 0

		wolfModel->draw(animShader);
		animShader->unbind();

		// Pop matrix stacks
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
	std::string resourceDir = "../resources";

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
	application->initParticleSystem(resourceDir);
	application->initGround();
	application->initWaterSurface(5.0f);
	application->initFoliage();
	application->initStars();
	application->initNature();
	application->initRocks();
	application->initRingRocks();

	if (application->wolfModel->loadModel(resourceDir + "/models/wolf.fbx")) {
		application->wolfModel->setScale(0.025f);

		// Start with animation 1 (hide|walk) and the appropriate speed
		application->wolfModel->setAnimation(1);
		application->wolfSpeed = 0.065;

		// Set materials for each mesh based on Blender materials
		// Mesh 0: Body (dark grey/brown fur)
		application->wolfModel->setMeshMaterial(0,
			glm::vec3(0.13f, 0.10f, 0.08f),  // Ambient
			glm::vec3(0.60f, 0.45f, 0.35f),  // Diffuse
			glm::vec3(0.3f, 0.3f, 0.3f),     // Specular
			16.0f);                          // Shininess

		// Mesh 1: Face/head (slightly lighter)
		application->wolfModel->setMeshMaterial(1,
			glm::vec3(0.14f, 0.11f, 0.09f),  // Ambient
			glm::vec3(0.65f, 0.50f, 0.40f),  // Diffuse
			glm::vec3(0.3f, 0.3f, 0.3f),     // Specular
			16.0f);                          // Shininess

		// Mesh 2: Eyes (dark)
		application->wolfModel->setMeshMaterial(2,
			glm::vec3(0.02f, 0.02f, 0.02f),  // Ambient
			glm::vec3(0.1f, 0.1f, 0.1f),     // Diffuse
			glm::vec3(0.5f, 0.5f, 0.5f),     // Specular
			32.0f);                          // Shininess

		// Mesh 3: Inner mouth (reddish)
		application->wolfModel->setMeshMaterial(3,
			glm::vec3(0.1f, 0.08f, 0.06f),   // Ambient
			glm::vec3(0.45f, 0.35f, 0.25f),  // Diffuse
			glm::vec3(0.3f, 0.3f, 0.3f),     // Specular
			16.0f);                          // Shininess

		// Mesh 4: Teeth (white)
		application->wolfModel->setMeshMaterial(4,
			glm::vec3(0.15f, 0.15f, 0.15f),  // Ambient
			glm::vec3(0.8f, 0.8f, 0.8f),     // Diffuse
			glm::vec3(0.5f, 0.5f, 0.5f),     // Specular
			32.0f);                          // Shininess

		// Mesh 5: Claws/paws (dark)
		application->wolfModel->setMeshMaterial(5,
			glm::vec3(0.15f, 0.05f, 0.05f),  // Ambient
			glm::vec3(0.7f, 0.3f, 0.3f),     // Diffuse
			glm::vec3(0.2f, 0.2f, 0.2f),     // Specular
			8.0f);                           // Shininess
	}
	else {
		std::cerr << "Failed to load wolf model!" << std::endl;
	}

	// compute sphere edge
	application->sphereEdge = application->sphereRadius * application->multiMeshes["sphereWTex"].scale * 2.0f - 2.0f;

	while (! glfwWindowShouldClose(windowManager->getHandle())) { // Loop until the user closes the window.
		application->render(); // Render scene
		glfwSwapBuffers(windowManager->getHandle()); // Swap front and back buffers
		glfwPollEvents(); // Poll for and process events
	}
	
	windowManager->shutdown(); // Quit program
	return 0;
}