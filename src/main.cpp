/**
 * Provided code for particle system simulator.
 * This code provides the mouse interface for clicking and dragging particles, and the
 * code to draw the system.  When the simulator is running system.advanceTime is called
 * to numerically integrate the system forward.
 * @author kry
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <sys/time.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Text.hpp"
#include "RecordVideo.hpp"
#include "UnifiedFramework.hpp"

using namespace std;

Scene cur_sim;
UnifiedFramework unified_sim;

bool run = false;
float stepsize = 1.0f / 60.0f; // 60 Hz screen refresh
int substeps = 1;

// parameters for interacting with particles
float maxDist = 150;
float minDist = 50;
float grabThresh = 10;

// for grabbing particles
Particle* p1 = NULL;
Particle* p2 = NULL;
float d1 = 0;
float d2 = 0;

int xdown = 0;
int ydown = 0;
int xcurrent = 0;
int ycurrent = 0;
bool mouseDown = false;
bool mouseInWindow = false;
bool grabbed = false;
bool wasPinned = false;
bool closeToParticlePairLine = false;
bool canEdit = true; // when simulation hasn't been run/stepped yet
bool gridShown = true;

string scene = "";
int testSceneID = 0;
int showKeyBindings = 0;
double avgTime = 0.0; 
double timeNum = 0.0;

// for openGL
GLFWwindow* window; // Main application window
string RES_DIR = ""; // Where data files live
shared_ptr<Program> progIM; // immediate mode

ImageRecorder imageRecorder;
bool recordFrames = false; // toggled by keyboard
bool recordFrame = false; // set to record when stepped in display

static void error_callback(int error, const char* description) {
	cerr << description << endl;
}
// Keyboard control pane
// TODO: Add some keys to control the parameters of the simulator
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {

	if ((action != GLFW_PRESS) && (action != GLFW_REPEAT)) return;
	if (key == GLFW_KEY_ESCAPE) {
		cout << avgTime/timeNum << endl;
		glfwSetWindowShouldClose(window, GL_TRUE);
	} else if (key == GLFW_KEY_SPACE) {
		// Space for running the simulator
		run = !run;
		canEdit = false;
	} else if (key == GLFW_KEY_O) {
		// o-O for gas_Open_Boundary scene
		cur_sim = GAS_OPEN;
        unified_sim.init(cur_sim);
	} else if (key == GLFW_KEY_C) {
		// c-C for Gas_Closed_Boundary scene
		cur_sim = GAS_CLOSED;
        unified_sim.init(cur_sim);
	} else if (key == GLFW_KEY_F) {
		// f-F for two-fluids scene
		cur_sim = TWO_FLUID;
		unified_sim.init(cur_sim);
	} else if (key == GLFW_KEY_Z) {
		// showing or unshowing grid lines
		if (gridShown) {
			unified_sim.unshowGrid();
			gridShown = false;
		} else {
			unified_sim.showGrid();
			gridShown = true;
		}
	} else if (key == GLFW_KEY_R && (mods & GLFW_MOD_SHIFT)) {
		// R for recording
		recordFrames = !recordFrames;
	} else if (key == GLFW_KEY_D) {
		showKeyBindings = (showKeyBindings+1)%3;
	} else if (key == GLFW_KEY_T) {
		// t-T for toggle gravity
		// TOdo : Add toggle gravity here
		// particleSystem.useGravity = !particleSystem.useGravity;
	} 
	
	// For managing substeps
	if (key == GLFW_KEY_UP) {
		substeps++;
	} else if (key == GLFW_KEY_DOWN) {
		substeps = max(substeps - 1, 1);
	}

	// [s-S] -> For making the simulator faster or slower
	if (key == GLFW_KEY_S) {
		if ((mods & GLFW_MOD_SHIFT)) {
			int oldsubsteps = substeps;
			substeps = max(1, substeps / 2);
			stepsize *= (float)substeps / (float)oldsubsteps;
		} else {
			int oldsubsteps = substeps;
			substeps = min(128, substeps * 2);
			stepsize *= (float)substeps / (float)oldsubsteps;
		}
	}

	float scale = (mods & GLFW_MOD_SHIFT) ? 1.01f : 1 / 1.01f;

	// [h,H] -> stepsize
	if (key == GLFW_KEY_H) {
		stepsize *= scale;
	} 
}

void mouse_pos_callback(GLFWwindow* window, double x, double y) {}

void mouse_callback(GLFWwindow* window, int button, int action, int mods) {}

static void init() {
	GLSL::checkVersion();

	// Check how many texture units are supported in the vertex shader
	int tmp;
	glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &tmp);
	cout << "GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS = " << tmp << endl;
	// Check how many uniforms are supported in the vertex shader
	glGetIntegerv(GL_MAX_VERTEX_UNIFORM_COMPONENTS, &tmp);
	cout << "GL_MAX_VERTEX_UNIFORM_COMPONENTS = " << tmp << endl;
	glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &tmp);
	cout << "GL_MAX_VERTEX_ATTRIBS = " << tmp << endl;

	// Set background color.
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// Enable z-buffer test.
	glEnable(GL_DEPTH_TEST);

	// Initialize the GLSL programs.
	progIM = make_shared<Program>();
	progIM->setVerbose(true);
	progIM->setShaderNames(RES_DIR + "simple_vert.glsl", RES_DIR + "simple_frag.glsl");
	progIM->init();
	progIM->addUniform("P");
	progIM->addUniform("MV");
	progIM->setVerbose(false);

	initTextRender(RES_DIR);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glDisable(GL_DEPTH_TEST);

	if (scene.empty()) {
		// The opening scene!
		// somthing here
	} else {
		// particleSystem.name = scene;
		// TODO: Add support for showing scene here
	}
	
	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

/** draws a line from the given point to the given particle */
void drawLineToParticle(double x, double y, Particle* p, double d) {
	if (p == NULL) return;
	if (d > maxDist) return;
	double col = d < minDist ? 1 : (maxDist - d) / (maxDist - minDist);
	glColor4d(1 - col, 0, col, 0.75f);
	glBegin(GL_LINES);
	glVertex2d(x, y);
	glVertex2d(p->p.x, p->p.y);
	glEnd();
}

void display() {
	struct timeval begin, end;
	if (run) {
		for (int i = 0; i < substeps; i++) {
			gettimeofday(&begin, 0);
			unified_sim.advanceTime(stepsize / substeps);
			gettimeofday(&end, 0);
			long seconds = end.tv_sec - begin.tv_sec;
			long microseconds = end.tv_usec - begin.tv_usec;
			double elapsed = seconds + microseconds*1e-6;
			avgTime += elapsed;
			timeNum += 1.0;
			//cout << elapsed << endl;
		}
		if ( recordFrames ) recordFrame = true;
	}

	unified_sim.display();
}

static void render() {
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	float aspect = width / (float)height;
	glViewport(0, 0, width, height);
	unified_sim.resize(width, height);

	// Clear framebuffer.
	//TODO: Set something fancy
	//glClearColor(0.2f, 0.2f, 0.1f, 0.0f);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-15, 15, -15, 15, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	display();

	stringstream ss;
	ss << "\n";
	if (showKeyBindings == 0) {
		ss << "Number of Particles: " << unified_sim.getNumParticles() << "\n";
		ss << "[f,F] Two-phase Fluid \n";
		ss << "[o,O] Gas Open Boundary \n";
		ss << "[c,C] Gas Closed Boundary \n";
		ss << "[z,Z] Draw Grid Lines \n";
		ss << "g Value = " << unified_sim.gravity[1] << "\n";
		ss << "[space] : start stop" << "\n";
	} 

	string text = ss.str();
	glm::mat4 projection = glm::ortho(0.0f, static_cast<float>(width), static_cast<float>(height), 0.0f);
	glm::mat4 modelview = glm::identity<glm::mat4>();

	RenderString(projection, modelview, 990, 50, 0.4, text);

	if (recordFrame) {
		recordFrame = false;
		imageRecorder.writeCurrentFrameToFile(window);
	}

	GLSL::checkError(GET_FILE_LINE);
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Please specify the resource directory." << endl;
		return 0;
	} else if (argc < 3) {
		cout << "Unspecified scene to run, using default value " << scene << endl;
	} else {
		scene = string(argv[2]);
	}
	RES_DIR = argv[1] + string("/");

	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if (!glfwInit()) {
		return -1;
	}
	// https://en.wikipedia.org/wiki/OpenGL
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	// glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(1280, 720, "COMP 559 W23 - Final Project - Mohanna Shahrad", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current
	glfwMakeContextCurrent(window);

	// Initialize GLEW.
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); 
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.

	glfwSwapInterval(1);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, mouse_pos_callback);

	// Initialize scene.
	init();

	// Loop until the user closes the window.
	while (!glfwWindowShouldClose(window)) {
		render();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// Quit program.
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}