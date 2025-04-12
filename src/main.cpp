#include <iostream>
#include <glad/glad.h>
#include "GLSL.h"
#include "Program.h"
#include "MatrixStack.h"
#include "WindowManager.h"
#include "Texture.h"
#define _USE_MATH_DEFINES
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Model.h"

using namespace std;
using namespace glm;

class Application : public EventCallbacks {
public:

	WindowManager* windowManager = nullptr;

	// shader programs
	std::shared_ptr<Program> prog, texProg, animProg;

	// ground data
	GLuint GrndBuffObj, GrndNorBuffObj, GIndxBuffObj;
	int g_GiboLen;
	GLuint GroundVertexArrayID;
	float groundSize = 20.0f;

	// wolf data
	Model* wolfModel;
	//Shape* wolfModel;
	float lastTime;
	float wolfSpeed = 5.0f;
	glm::vec3 wolfPosition = glm::vec3(0.0f, 0.0f, 0.0f);
	float wolfRotation = 0.0f;

	// plant data
	//Shape* plantModel;

	// camera variables
	float phi = 0.0f;                               // Pitch angle in radians
	float theta = glm::pi<float>();                 // Yaw angle in radians
	float radius = 5.0f;                            // Distance from camera to look-at point
	glm::vec3 eye = glm::vec3(0.0f, 3.0f, 5.0f);    // Camera position
	glm::vec3 lookAt = glm::vec3(0.0f, 0.0f, 0.0f); // Point to look at
	double lastX, lastY;
	bool mouseInitialized = false;
	float mouseSpeed = 0.005f;                      // Mouse sensitivity

	// input state
	bool moveForward = false;
	bool moveBackward = false;
	bool moveLeft = false;
	bool moveRight = false;

	void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
			glfwSetWindowShouldClose(window, GL_TRUE);
		}

		// WASD movement keys
		if (key == GLFW_KEY_W) {
			moveForward = (action == GLFW_PRESS || action == GLFW_REPEAT);
		}
		if (key == GLFW_KEY_S) {
			moveBackward = (action == GLFW_PRESS || action == GLFW_REPEAT);
		}
		if (key == GLFW_KEY_A) {
			moveLeft = (action == GLFW_PRESS || action == GLFW_REPEAT);
		}
		if (key == GLFW_KEY_D) {
			moveRight = (action == GLFW_PRESS || action == GLFW_REPEAT);
		}

		// Animation debug
		if (key == GLFW_KEY_N && action == GLFW_PRESS) {
			int nextAnim = (wolfModel->getCurrentAnimation() + 1) % wolfModel->getAnimationCount();
			wolfModel->setAnimation(nextAnim);
		}
		if (key == GLFW_KEY_P && action == GLFW_PRESS) {
			int prevAnim = wolfModel->getCurrentAnimation() - 1;
			if (prevAnim < 0) prevAnim = wolfModel->getAnimationCount() - 1;
			wolfModel->setAnimation(prevAnim);
		}
		if (key == GLFW_KEY_0 && action == GLFW_PRESS && mods == GLFW_MOD_CONTROL) {
			std::cout << "Available animations:" << std::endl;
			for (int i = 0; i < wolfModel->getAnimationCount(); i++) {
				std::cout << "  " << i << ": " << wolfModel->getAnimationName(i) << std::endl;
			}
		}
	}

	void mouseCallback(GLFWwindow* window, int button, int action, int mods) override {
		// Implement if needed for mouse button clicks
	}

	void scrollCallback(GLFWwindow* window, double deltaX, double deltaY) override {
		// Scroll callback is now empty as we don't want the mouse to affect the camera
	}

	void mouseMoveCallback(GLFWwindow* window, double xpos, double ypos) {
		if (!mouseInitialized) {
			lastX = xpos;
			lastY = ypos;
			mouseInitialized = true;
			return;
		}

		float deltaX = xpos - lastX;
		lastX = xpos;
		lastY = ypos;

		// Mouse only affects wolf rotation now
		wolfRotation -= deltaX * mouseSpeed;  // Negative for intuitive control

		// Keep the rotation in the range [0, 2π]
		if (wolfRotation > 2 * glm::pi<float>()) wolfRotation -= 2 * glm::pi<float>();
		if (wolfRotation < 0) wolfRotation += 2 * glm::pi<float>();

		updateCamera();
	}

	void updateCamera() {
		// Set fixed camera height and distance
		float cameraHeight = 3.0f;   // Height above the wolf
		float cameraDistance = 5.0f; // Distance behind the wolf

		// Calculate forward direction of the wolf
		glm::vec3 wolfForward = glm::vec3(
			sin(wolfRotation),  // x component
			0.0f,              // y component
			cos(wolfRotation)   // z component
		);

		// Position camera behind and above the wolf
		eye = wolfPosition - (wolfForward * cameraDistance) + glm::vec3(0.0f, cameraHeight, 0.0f);

		// Look slightly ahead of the wolf's position
		lookAt = wolfPosition + (wolfForward * 3.0f);
	}

	void resizeCallback(GLFWwindow* window, int width, int height) {
		glViewport(0, 0, width, height);
	}

	void init(const std::string& resourceDirectory) {
		GLSL::checkVersion();
		glClearColor(0.5f, 0.7f, 1.0f, 1.0f);      // Set background color to sky blue
		glEnable(GL_DEPTH_TEST);                    // Enable depth testing

		// Initialize basic shader
		prog = make_shared<Program>();
		prog->setVerbose(true);
		prog->setShaderNames(resourceDirectory + "/simple_vert.glsl", resourceDirectory + "/simple_frag.glsl");
		prog->init();
		prog->addUniform("P");
		prog->addUniform("V");
		prog->addUniform("M");
		prog->addUniform("MatAmb");
		prog->addUniform("MatDif");
		prog->addUniform("MatSpec");
		prog->addUniform("MatShine");
		prog->addUniform("lightPos");
		prog->addAttribute("vertPos");
		prog->addAttribute("vertNor");

		// Initialize texture shader
		texProg = make_shared<Program>();
		texProg->setVerbose(true);
		texProg->setShaderNames(resourceDirectory + "/tex_vert.glsl", resourceDirectory + "/tex_frag.glsl");
		texProg->init();
		texProg->addUniform("P");
		texProg->addUniform("V");
		texProg->addUniform("M");
		texProg->addUniform("Texture");
		texProg->addUniform("lightPos");
		texProg->addAttribute("vertPos");
		texProg->addAttribute("vertNor");
		texProg->addAttribute("vertTex");

		// Initialize animation shader
		animProg = make_shared<Program>();
		animProg->setVerbose(true);
		animProg->setShaderNames(resourceDirectory + "/anim_vert.glsl", resourceDirectory + "/anim_frag.glsl");
		animProg->init();
		animProg->addUniform("P");
		animProg->addUniform("V");
		animProg->addUniform("M");

		// Material uniforms
		animProg->addUniform("material.ambient");
		animProg->addUniform("material.diffuse");
		animProg->addUniform("material.specular");
		animProg->addUniform("material.emissive");
		animProg->addUniform("material.shininess");
		animProg->addUniform("material.opacity");
		animProg->addUniform("material.hasTexture");


		animProg->addUniform("lightPos");

		// Add uniforms for bone transforms
		for (int i = 0; i < 100; i++) {
			animProg->addUniform("boneTransforms[" + std::to_string(i) + "]");
		}

		// Add uniforms for textures
		animProg->addUniform("texture_diffuse1");
		animProg->addUniform("texture_specular1");
		animProg->addUniform("texture_normal1");
		animProg->addUniform("texture_height1");
		animProg->addUniform("texture_emissive1");

		// Add attributes
		animProg->addAttribute("vertPos");
		animProg->addAttribute("vertNor");
		animProg->addAttribute("vertTex");
		animProg->addAttribute("vertTan");
		animProg->addAttribute("vertBitan");
		animProg->addAttribute("boneIDs");
		animProg->addAttribute("weights");


		// Create the wolf model
		wolfModel = new Model();
		//wolfModel = new Shape();
		lastTime = glfwGetTime();

		//plantModel = new Shape();
	}

	// Initialize the ground plane
	void initGround() {
		float g_groundSize = groundSize;
		float g_groundY = 0.0f;
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
		unsigned short idx[] = { 0, 1, 2, 0, 2, 3 };

		// Generate the ground VAO
		glGenVertexArrays(1, &GroundVertexArrayID);
		glBindVertexArray(GroundVertexArrayID);

		g_GiboLen = 6;
		glGenBuffers(1, &GrndBuffObj);
		glBindBuffer(GL_ARRAY_BUFFER, GrndBuffObj);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GrndPos), GrndPos, GL_STATIC_DRAW);

		glGenBuffers(1, &GrndNorBuffObj);
		glBindBuffer(GL_ARRAY_BUFFER, GrndNorBuffObj);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GrndNorm), GrndNorm, GL_STATIC_DRAW);

		glGenBuffers(1, &GIndxBuffObj);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GIndxBuffObj);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(idx), idx, GL_STATIC_DRAW);
	}

	// Set material properties
	void SetMaterial(shared_ptr<Program> curS, int i) {
		switch (i) {
		case 0: // Ground material (green grass color)
			glUniform3f(curS->getUniform("MatAmb"), 0.1, 0.2, 0.1);
			glUniform3f(curS->getUniform("MatDif"), 0.2, 0.6, 0.2);
			glUniform3f(curS->getUniform("MatSpec"), 0.1, 0.1, 0.1);
			glUniform1f(curS->getUniform("MatShine"), 5.0);
			break;
		case 1: // Wolf material
			glUniform3f(curS->getUniform("MatAmb"), 0.13f, 0.1f, 0.08f);
			glUniform3f(curS->getUniform("MatDif"), 0.6f, 0.45f, 0.35f);
			glUniform3f(curS->getUniform("MatSpec"), 0.3f, 0.3f, 0.3f);
			glUniform1f(curS->getUniform("MatShine"), 16.0f);
			break;
		}
	}

	// Draw the ground plane using matrix stack
	void drawGround(shared_ptr<Program> curS, std::shared_ptr<MatrixStack> Model) {
		curS->bind();
		glBindVertexArray(GroundVertexArrayID);

		// Set material for ground
		SetMaterial(curS, 0);

		// Use the matrix stack for the model matrix
		Model->pushMatrix();
		Model->loadIdentity();
		Model->translate(vec3(0, 0, 0));
		glUniformMatrix4fv(curS->getUniform("M"), 1, GL_FALSE, value_ptr(Model->topMatrix()));

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, GrndBuffObj);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, GrndNorBuffObj);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);

		// Draw
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GIndxBuffObj);
		glDrawElements(GL_TRIANGLES, g_GiboLen, GL_UNSIGNED_SHORT, 0);

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		Model->popMatrix();
		curS->unbind();
	}

	// Update the wolf based on user input
	void updateWolf(float dt) {
		// Initialize movement direction
		float forwardAmount = 0.0f;
		float rightAmount = 0.0f;

		// Standard WASD mapping
		if (moveForward)
			forwardAmount += 1.0f;  // W - forward
		if (moveBackward)
			forwardAmount -= 1.0f;  // S - backward
		if (moveLeft)
			rightAmount -= 1.0f;    // A - left
		if (moveRight)
			rightAmount += 1.0f;    // D - right

		// Check if there's any movement input
		bool isMoving = (forwardAmount != 0.0f || rightAmount != 0.0f);

		if (isMoving) {
			// Calculate forward and right vectors based on wolf's rotation
			glm::vec3 forward = glm::vec3(
				sin(wolfRotation),  // x component
				0.0f,              // y component
				cos(wolfRotation)   // z component
			);

			// Right vector is perpendicular to forward
			glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), forward));

			// Use the existing wolfSpeed variable
			float moveAmount = wolfSpeed * dt;

			// Move the wolf based on input
			wolfPosition += forward * moveAmount * forwardAmount; // Forward/backward
			wolfPosition += right * moveAmount * -rightAmount;     // Left/right

			// Set walk animation when moving
			if (wolfModel->getCurrentAnimation() != 1) {
				wolfModel->setAnimation(1);
			}
		}
		else {
			// Set idle animation when not moving
			if (wolfModel->getCurrentAnimation() != 0) {
				wolfModel->setAnimation(0);
			}
		}

		// Update wolf model position and rotation
		wolfModel->setPosition(wolfPosition);

		// Apply a -180 degree rotation when drawing the wolf to align it with the movement direction
		wolfModel->setRotation(wolfRotation - glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));
	}

	// Draw the wolf using matrix stack
	void drawWolf(shared_ptr<Program> curS, std::shared_ptr<MatrixStack> Model) {

		// Use matrix stack for wolf model
		Model->pushMatrix();
		Model->loadIdentity();
		curS->bind();

		// Set base material properties (these will be overridden by per-mesh materials)
		glUniform3f(curS->getUniform("material.ambient"), 0.1f, 0.1f, 0.1f);
		glUniform3f(curS->getUniform("material.diffuse"), 0.8f, 0.8f, 0.8f);
		glUniform3f(curS->getUniform("material.specular"), 0.3f, 0.3f, 0.3f);
		glUniform3f(curS->getUniform("material.emissive"), 0.0f, 0.0f, 0.0f);
		glUniform1f(curS->getUniform("material.shininess"), 32.0f);
		glUniform1f(curS->getUniform("material.opacity"), 1.0f);

		// Set texture uniforms

		// The Model class already handles position and rotation internally,
		// but we can use the matrix stack to add additional transformations if needed
		// For now, we'll just use an identity matrix
		glUniformMatrix4fv(curS->getUniform("M"), 1, GL_FALSE, value_ptr(Model->topMatrix()));

		// Draw the wolf model
		wolfModel->draw(curS);

		Model->popMatrix();
		curS->unbind();
	}

	void render() {
		int width, height;
		glfwGetFramebufferSize(windowManager->getHandle(), &width, &height);
		glViewport(0, 0, width, height);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_PROGRAM_POINT_SIZE);

		float aspect = width / (float)height;

		// Calculate delta time for smooth animation and movement
		float currentTime = glfwGetTime();
		float dt = currentTime - lastTime;
		lastTime = currentTime;

		// Update wolf position and animation based on user input
		updateWolf(dt);

		// Update camera position to follow the wolf
		updateCamera();

		// Update wolf animation
		wolfModel->update(dt);

		// Create matrix stacks
		auto Projection = make_shared<MatrixStack>();
		auto View = make_shared<MatrixStack>();
		auto Model = make_shared<MatrixStack>();

		Projection->pushMatrix();
		Projection->perspective(45.0f, aspect, 0.01f, 100.0f);

		View->pushMatrix();
		View->loadIdentity();

		// Set view matrix with updated camera
		glm::mat4 viewMatrix = glm::lookAt(eye, lookAt, glm::vec3(0.0f, 1.0f, 0.0f));
		View->multMatrix(viewMatrix);

		// Lighting position (above the scene)
		glm::vec3 lightPos(0.0f, 10.0f, 0.0f);

		// Draw the ground with matrix stack
		prog->bind();
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
		glUniformMatrix4fv(prog->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));
		glUniform3f(prog->getUniform("lightPos"), lightPos.x, lightPos.y, lightPos.z);
		drawGround(prog, Model);
		prog->unbind();

		// Draw the wolf with matrix stack
		animProg->bind();
		glUniformMatrix4fv(animProg->getUniform("P"), 1, GL_FALSE, value_ptr(Projection->topMatrix()));
		glUniformMatrix4fv(animProg->getUniform("V"), 1, GL_FALSE, value_ptr(View->topMatrix()));
		glUniform3f(animProg->getUniform("lightPos"), lightPos.x, lightPos.y, lightPos.z);
		drawWolf(animProg, Model);
		animProg->unbind();

		prog->bind();

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

int main(int argc, char* argv[]) {
	std::string resourceDir = "../resources";

	if (argc >= 2) {
		resourceDir = argv[1];
	}

	Application* application = new Application();
	WindowManager* windowManager = new WindowManager();
	windowManager->init(1280, 720);
	windowManager->setEventCallbacks(application);
	application->windowManager = windowManager;

	glfwSetWindowUserPointer(windowManager->getHandle(), application);
	glfwSetCursorPosCallback(windowManager->getHandle(), mouseMoveCallbackWrapper);
	glfwSetScrollCallback(windowManager->getHandle(), scrollCallbackWrapper);

	application->init(resourceDir);
	application->initGround();

	// Load the wolf model
	if (application->wolfModel->loadModel(resourceDir + "/models/wolf.fbx")) { // Changed from loadModel to loadMesh
		application->wolfModel->setScale(glm::vec3(0.025f));
		application->wolfModel->setAnimation(1);
	}
	else {
		std::cerr << "Failed to load wolf model!" << std::endl;
	}

	while (!glfwWindowShouldClose(windowManager->getHandle())) {
		application->render();
		glfwSwapBuffers(windowManager->getHandle());
		glfwPollEvents();
	}

	windowManager->shutdown();
	return 0;
}