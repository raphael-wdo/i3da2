/*
	NAME: I3D ASSINGMENT 2 2020 SEM 1
	AUTHOR: RAPHAEL DOH ONG WONG (S3735236)
	DESC: COSC1186/1187 Interactive 3D Graphics and Animation Assignment 2 - 3D Island Survival Game
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#if _WIN32
#   include <Windows.h>
#endif
#if __APPLE__
#   include <OpenGL/gl.h>
#   include <OpenGL/glu.h>
#   include <GLUT/glut.h>
#else
#   include <GL/gl.h>
#   include <GL/glu.h>
#   include <GL/glut.h>
#endif

#if 0
// Program uses the Simple OpenGL Image Library for loading textures: http://www.lonesock.net/soil.html
#include <SOIL.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

using namespace std;

GLuint texture;

typedef struct { float x, y; } vec2f;
vec2f mouseZoom = { 0,0 };

typedef struct { float x, y, x0, y0, zoom; } camera;
camera c = { 330, 35, 0, 0, -0.5 };

int mouseMode = 0;
GLboolean wireframe = FALSE;
GLboolean normalview = FALSE;
GLboolean tangentview = FALSE;
GLboolean lightview = TRUE;
GLboolean textureview = TRUE;
GLboolean axes = TRUE;
GLboolean wPressed, aPressed, sPressed, dPressed, leftKeyPressed, rightKeyPressed = FALSE;
GLboolean gameOver = FALSE;

//perspective vars
float fovy = 75;
float aspect = 0;
float zNear = 0.1;
float zFar = 100;

//number of rows
int n = 64;

//wave vars
typedef struct { float A, k, kz, w; } sinewave;
sinewave sw[] = {
	{ 0.02, 8 * M_PI, 4 * M_PI, 0.25 * M_PI },
	{ 0.02, 4 * M_PI, 8 * M_PI, 0.5 * M_PI }
};
float wavelength = 0.5;
float waveheight = 0.06;

//time vars
typedef struct { float t, lastT, dt; } time;
time waveTime = { 0.0 };

//light var
float pos[] = { 0.9,10,0.9,0 };
float white[] = { 1.0,1.0,1.0,0 };
float blue[] = { 0.12,0.12,0.4,0 };
float green[] = { 0.0,0.4,0.0,0 };
float red[] = { 0.4,0.0,0.0,0 };
float gray[] = { 0.2,0.2,0.2,0 };
float shine[] = { 128,0,0,0 };

//boat vars
typedef struct { float x, z; GLboolean hit; } boat;
boat enemies[128];
int boatCount = 0;
float boatSpawnTimer = 3;
int score = 0;

//cannon vars
vec2f cannon_angle = {60.0, 0.0};
typedef struct { float speed, angle; } vec2fPolar;
vec2fPolar initVel = { 0.5, 60.0 };
vec2fPolar boatVel = { 0.6, 60.0 };
typedef struct { float x, y, z; } vec3f;
typedef struct { vec3f r0, v0, r, v; float angle, phi, launchedTime, cnnPower; GLboolean sunk; } state;
state cannonballs[10000];
int cnnbllCount = 0;
float shootTimer = 2;

//island vars
float island_r = 0.25;
int health = 100;

//performance info var
typedef struct {
	GLboolean debug;
	GLboolean go;
	float startTime;
	GLboolean OSD;
	int frames;
	float frameRate;
	float frameRateInterval;
	float lastFrameRateT;
} global_t;

global_t global = { TRUE, TRUE, 0.0, TRUE, 0, 0.0, 0.2, 0.0 };
float lastT = -1.0;

//skybox
typedef struct { GLuint xn, xp, yn, yp, zn, zp; } skyboxTextures;
skyboxTextures skyboxTexture;

const float g = -0.25;
const int milli = 1000;

float rad2deg(float rad)
{
	return rad * 180 / M_PI;
}
float deg2rad(float deg)
{
	return M_PI / 180 * deg;
}

// load a texture from file using the stb_image library
uint32_t loadTexture(const char* filename) {
	int width, height, components;
	unsigned char* data = stbi_load(filename, &width, &height, &components, STBI_rgb);

	glPushAttrib(GL_TEXTURE_BIT);

	unsigned int id;
	glGenTextures(1, &id);
	glBindTexture(GL_TEXTURE_2D, id);
	//      glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

	glPopAttrib();

	printf("texture loaded\n");
	return id;
}


#if 0
GLuint loadTexture(const char* filename)
{
	GLuint tex = SOIL_load_OGL_texture(filename, SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y);
	if (!tex)
		return 0;

	glBindTexture(GL_TEXTURE_2D, tex);
#if 0
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
#endif
	glBindTexture(GL_TEXTURE_2D, 0);

	return tex;
}
#endif

void drawAxes(float len)
{
#if 0
	glPushAttrib(GL_CURRENT_BIT);
#endif
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(-len, 0.0, 0.0);
	glVertex3f(len, 0.0, 0.0);
	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, -len, 0.0);
	glVertex3f(0.0, len, 0.0);
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, -len);
	glVertex3f(0.0, 0.0, len);
	glEnd();
#if 0
	glPopAttrib();
#endif
}

void cameraFunc(void)
{
	//inverse values of angles
	glTranslatef(0.2, -0.15, c.zoom);
	glRotatef(c.y, 1, 0, 0);
	glRotatef(c.x -45, 0, 1, 0);

}

void inverseCamFunc()
{
	glTranslatef(-0.2, 0.15, c.zoom);
	glRotatef(360 - c.y, 1, 0, 0);
	glRotatef(360 - c.x - 45, 0, 1, 0);
}

void calcSineWave(sinewave sw, float x, float z, float t, float* y, GLboolean dvs, float* dydx, float* dydz)
{
	float angle = sw.k * x + sw.kz * z + sw.w * t;
	*y = sw.A * sinf(angle);
	if (dvs) {
		*dydx = sw.A * sw.k * cosf(angle);
		*dydz = sw.A * sw.kz * cosf(angle);
	}
}

void drawNormal(float x, float z, float t)
{
	float dydx = 0;
	float dydz = 0;

	//find height of sine wave
	float y = 0;
	float sum = 0;
	for (int wave = 0; wave < 2; wave++) {
		calcSineWave(sw[wave], x, z, t, &y, TRUE, &dydx, &dydz);
		sum += y;
	}
	y = sum;

	//normal and tangent

	if (normalview) {
		glColor3f(1, 1, 0);
		glVertex3f(x, y, z);
		glVertex3f(x - 0.1 * dydx, y + 0.1, z - 0.1 * dydz);
	}
	if (tangentview) {
		glColor3f(0.0, 1.0, 1.0);
		glVertex3f(x, y, z);
		glVertex3f(x + 0.1, y + 0.1 * dydx, z);
		glVertex3f(x, y, z);
		glVertex3f(x, y + 0.1 * dydz, z + 0.1);
	}

	glColor3f(1.0, 1.0, 1.0);

}

void setLight()
{
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, white);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
	glLightfv(GL_LIGHT0, GL_SPECULAR, white);
}

void setMaterial(float* color) 
{
	setLight();
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
}

void setVertex(float x, float y, float z)
{
	glNormal3f(x,y + 0.01 ,z);
	glVertex3f(x, y, z);
}

void drawGrid()
{
	setMaterial(blue);
	glColor3f(0.0, 0.0, 1.0);

	float xStep = 2.0 / n;
	float zStep = 2.0 / n;
	for (int j = 0; j < n; j++) {
		glBegin(GL_TRIANGLE_STRIP);
		float z = -1.0 + j * zStep;
		float t = waveTime.t;

		//draw Grid
		for (int i = 0; i <= n; i++) {
			float x = -1.0 + i * xStep;

			//find height of sine wave
			float y = 0;
			float sum = 0;
			float dydx = 0;
			float dydz = 0;

			for (int wave = 0; wave < 2; wave++){
				calcSineWave(sw[wave], x, z, t, &y, FALSE, &dydx, &dydz);
				sum += y;
			}
			y = sum;
			sum = 0;	//MUST RESET SUM OR ELSE GAP WILL APPEAR DUE TO Y BEING DOUBLED
			setVertex(x, y, z);

			float z1 = -1.0 + ((j + 1) * zStep);
			for (int wave = 0; wave < 2; wave++) {
				calcSineWave(sw[wave], x, z1, t, &y, FALSE, &dydx, &dydz);
				sum += y;
			}
			y = sum;
			setVertex(x, y, z1);
		}
		glEnd();

		//draw Normals
		if (normalview || tangentview) {
			glBegin(GL_LINES);
			for (int i = 0; i <= n; i++) {
				float x = -1.0 + i * xStep;
				drawNormal(x, z, t);
				float z1 = -1.0 + ((j + 1) * zStep);
				drawNormal(x, z1, t);
			}
			glEnd();
		}
	}

	glColor3f(1.0, 1.0, 1.0);
}

void drawGun(float angle, float r, float h)
{
	const int slices = n;
	float theta;
	float x1, y1, z1, y2;
	glRotatef(angle-90, 1, 0, 0);
	glBegin(GL_QUAD_STRIP);
	for (int i = 0; i <= slices; i++) {
		theta = i / (float)slices * 2.0 * M_PI;
		x1 = r * sinf(theta);
		y1 = 0.0;
		z1 = r * cosf(theta);
		y2 = h;
		setVertex(x1, y1, z1);
		setVertex(x1, y2, z1);
	}
	glEnd();

	glRotatef(360 - angle + 90, 1, 0, 0);
}



void displayVelocity(vec2f cannon_angle, float hull_angle)
{
	glTranslatef(0.0, 0.1, 0.0);
	glRotatef(cannon_angle.y + 180, 0, 1, 0);
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	initVel.angle = cannon_angle.x;

	float x = initVel.speed * cosf(deg2rad(initVel.angle + hull_angle));
	float y = initVel.speed * sinf(deg2rad(initVel.angle + hull_angle));

	glVertex3f(x, y, 0.0);
	glEnd();
	glRotatef(360 - cannon_angle.y + 180, 0, 1, 0);
	glTranslatef(0.0, -0.1, 0.0);

}


void drawTrajectoryAnalytical(vec2f cannon_angle, float hull_angle, float time)
{
	setMaterial(white);
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 1.0);

	float y = 0;

	// define projectile location across the cannon barrel
	vec2f r, cannon_r;
	vec3f rr;

	//calculation formulas
	float phi = deg2rad(cannon_angle.y);
	float x_dir = cosf(deg2rad(cannon_angle.x + hull_angle));
	float y_dir = sinf(deg2rad(cannon_angle.x + hull_angle));
	cannon_r.x = 0.15 * x_dir;
	cannon_r.y = 0.15 * y_dir + 0.1;

	for (float t = time; y >= 0; t += 0.1)
	{
		//cant use glut as the values are too small and cause many unncessary loops

		// x = vector speed in direction x * time + origin position
		r.x = initVel.speed * x_dir * t + cannon_r.x;

		// y = -1/2 * gravity * time^2 + vector speed in direction y * time + origin position
		r.y = (-0.5 * -g * pow(t, 2)) + (initVel.speed * y_dir * t) + cannon_r.y;
		y = r.y;

		rr.x = r.x * cosf(phi);
		rr.z = r.x * -sinf(phi);
		rr.y = r.y;

		glVertex3f(rr.x, rr.y, rr.z);
		//printf("x: %.2f\ty: %.2f\ttime: %.2f\n", x, y, time);

	}
	glEnd();
}

void drawBallTrajectory(state projectile)
{
	setMaterial(red);
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 1.0);

	float y = 0;

	// define projectile location across the cannon barrel
	vec2f r, cannon_r;
	vec3f rr;

	//calculation formulas
	float phi = deg2rad(projectile.phi);
	float x_dir = cosf(deg2rad(projectile.angle));
	float y_dir = sinf(deg2rad(projectile.angle));
	cannon_r.x = 0.15 * x_dir;
	cannon_r.y = 0.15 * y_dir + 0.1;

	for (float time = waveTime.t - projectile.launchedTime; y >= 0; time += 0.1)
	{
		//cant use glut as the values are too small and cause many unncessary loops

		// x = vector speed in direction x * time + origin position
		r.x = projectile.cnnPower * x_dir * time + cannon_r.x;

		// y = -1/2 * gravity * time^2 + vector speed in direction y * time + origin position
		r.y = (-0.5 * -g * pow(time, 2)) + (projectile.cnnPower * y_dir * time) + cannon_r.y;
		y = r.y;

		rr.x = r.x * cosf(phi) + projectile.r0.x;
		rr.z = r.x * -sinf(phi) + projectile.r0.z;
		rr.y = r.y;

		glVertex3f(rr.x, rr.y, rr.z);
		//printf("x: %.2f\ty: %.2f\ttime: %.2f\n", x, y, time);

	}
	glEnd();
}


void updateProjectileStateAnalytical(float t, state projectile, int cnnbllNo, float hull_angle, float cannonPower)
{
	// define projectile location across the cannon barrel
	vec2f r, cannon_r;

	//calculation formulas
	float phi = deg2rad(projectile.phi);
	float x_dir = cosf(deg2rad(projectile.angle + hull_angle));
	float y_dir = sinf(deg2rad(projectile.angle + hull_angle));
	cannon_r.x = 0.15 * x_dir;
	cannon_r.y = 0.15 * y_dir + 0.1;

	// x = vector speed in direction x * time + origin position
	r.x = cannonPower * x_dir * t + cannon_r.x;

	// y = -1/2 * gravity * time^2 + vector speed in direction y * time + origin position
	r.y = (-0.5 * -g * pow(t, 2)) + (cannonPower * y_dir * t) + cannon_r.y;

	cannonballs[cnnbllNo].r.x = projectile.r.x = r.x * cosf(phi) + projectile.r0.x;
	cannonballs[cnnbllNo].r.z = projectile.r.z = r.x * -sinf(phi) + projectile.r0.z;
	cannonballs[cnnbllNo].r.y = projectile.r.y = r.y;

	//printf("inside x: %.2f\n", projectile.r.x);

}

void drawBall()
{
	for (int i = 0; i < cnnbllCount; i++)
	{
		state cnnbll = cannonballs[i];

		if (!cnnbll.sunk && cnnbll.r.y > -0.2) {
			glTranslatef(cnnbll.r.x, cnnbll.r.y, cnnbll.r.z);
			glutSolidSphere(0.01, 16, 16);
			glTranslatef(-cnnbll.r.x, -cnnbll.r.y, -cnnbll.r.z);
		}

		drawBallTrajectory(cnnbll);
	}

}

void drawSkyBox()
{
	setMaterial(white);
	glColor3f(1.0, 1.0, 1.0);

	glEnable(GL_TEXTURE_2D);

	int numOfSides = 6;
	int numOfVertex = 4;
	GLuint textureFactory[numOfSides] = {skyboxTexture.xn, skyboxTexture.xp, skyboxTexture.yn, skyboxTexture.yp, skyboxTexture.zn, skyboxTexture.zp, };
	vec3f vertexFactory[] = {
		{-1,1,-1},	//xn
		{-1,1,1},	
		{-1,-1,1},	
		{-1,-1,-1},	
	
		{1,1,1},	//xp
		{1,1,-1},
		{1,-1,-1},
		{1,-1,1},

		{-1,-1,-1},	//yn
		{1,-1,-1},
		{1,-1,1},
		{-1,-1,1},

		{-1,1,-1},	//yp
		{1,1,-1},
		{1,1,1},
		{-1,1,1},

		{1,1,-1},	//zn
		{-1,1,-1},
		{-1,-1,-1},
		{1,-1,-1},
		
		{-1,1,1},	//zp
		{1,1,1},
		{1,-1,1},
		{-1,-1,1},
	};
	vec2f texCoordFactory[] = {
		{0,0},
		{1,0},
		{1,1},
		{0,1}
	};
	
	for (int i = 0; i < numOfSides; i++) {
		glBindTexture(GL_TEXTURE_2D, textureFactory[i]);

		glBegin(GL_QUADS);
		for (int j = 0; j < numOfVertex; j++) {
			vec2f texCoord = texCoordFactory[j];
			vec3f vertex = vertexFactory[i * numOfVertex + j];
			glTexCoord2f(texCoord.x, texCoord.y);
			glVertex3f(vertex.x, vertex.y, vertex.z);
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);

	glColor3f(1.0, 1.0, 1.0);
}

void drawIsland()
{
	
	if (textureview) {
		setMaterial(white);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, texture);
	}
	else {
		setMaterial(green);
		glDisable(GL_TEXTURE_2D);
	}

	glColor3f(0.0, 1.0, 0.0);
	

	const int slices = n, stacks = n;
	float r = 0.3;	//island size
	float theta, phi;
	float x1, y1, z1, x2, y2, z2;
	float step_phi = M_PI / stacks;

	glTranslatef(0, -0.21, 0);

	for (int j = 0; j <= stacks / 2; j++) {
		phi = j / (float)stacks * M_PI;
		glBegin(GL_QUAD_STRIP);
		for (int i = 0; i <= slices; i++) {
			theta = i / (float)slices * 2.0 * M_PI;
			x1 = r * sinf(phi) * cosf(theta);
			y1 = r * cosf(phi); 
			z1 = r * sinf(phi) * sinf(theta);
			x2 = r * sinf(phi + step_phi) * cosf(theta);
			y2 = r * cosf(phi + step_phi);
			z2 = r * sinf(phi + step_phi) * sinf(theta);
			glTexCoord2f(x1, z1);
			setVertex(x1, y1, z1);
			glTexCoord2f(x2, z2);
			setVertex(x2, y2, z2);
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);

	glTranslatef(0, 0.2, 0);
	r = 0.125;
	setMaterial(gray);

	//cannon base
	glBegin(GL_QUAD_STRIP);
	for (int i = 0; i <= slices; i++) {
		theta = i / (float)slices * 2.0 * M_PI;
		x1 = r * sinf(theta);
		y1 = 0.0;
		z1 = r * cosf(theta);
		y2 = 0.1;
		setVertex(x1, y1, z1);
		setVertex(x1, y2, z1);
	}
	for (int i = 0; i <= slices; i++) {
		theta = i / (float)slices * 2.0 * M_PI;
		x1 = r * sinf(theta);
		z1 = r * cosf(theta);
		y2 = 0.1;
		setVertex(0, y2, 0);
		setVertex(x1, y2, z1);
	}
	glEnd();
	
	//cannon body
	r = 0.075;
	glTranslatef(0, 0.1, 0);
	for (int j = 0; j <= stacks / 2; j++) {
		phi = j / (float)stacks * M_PI;
		glBegin(GL_QUAD_STRIP);
		for (int i = 0; i <= slices; i++) {
			theta = i / (float)slices * 2.0 * M_PI;
			x1 = r * sinf(phi) * cosf(theta);
			y1 = r * cosf(phi);
			z1 = r * sinf(phi) * sinf(theta);
			x2 = r * sinf(phi + step_phi) * cosf(theta);
			y2 = r * cosf(phi + step_phi);
			z2 = r * sinf(phi + step_phi) * sinf(theta);
			setVertex(x1, y1, z1);
			setVertex(x2, y2, z2);
		}
		glEnd();
	}

	//cannon gun
	glRotatef(cannon_angle.y + 90 , 0, 1, 0);
	drawGun(cannon_angle.x,0.01,0.15);
	
	glRotatef(360-cannon_angle.y + 90, 0, 1, 0);
	glTranslatef(0, -0.09, 0); //reset pos

	glColor3f(1.0, 1.0, 1.0);

}

void drawBoat(float x, float z)
{
	setMaterial(red);
	glColor3f(1.0, 0.0, 0.0);

	//find height of sine wave
	float y = 0;
	float sum = 0;
	float dydx = 0;
	float dydz = 0;

	for (int wave = 0; wave < 2; wave++) {
		calcSineWave(sw[wave], x, z, waveTime.t, &y, TRUE, &dydx, &dydz);
		sum += y;
	}
	y = sum;

	glTranslatef(x, y, z);

	//find rotation
	float heading = rad2deg(atan(x / z));
	if (z >= 0) heading -= 180;
	glRotatef(heading, 0, 1, 0);
	glRotatef(rad2deg(dydx), 0, 0, 1);
	glRotatef(rad2deg(dydz), 1, 0, 0);

	int num_of_vectors = 18;
	float vecx[] = { -0.05,-0.05,-0.05,-0.05,0.05,0.05,0.05,0.05,-0.05,-0.05,-0.05,-0.05,0.05,0.05,0.05,0.05,-0.05,-0.05 };
	float vecy[] = { 0.05,0.05,0.0,0.0,0.0,0.0,0.05,0.05,0.05,0.05,0.05,0.0,0.05,0.0,0.0,0.05,0.0,0.05 };
	float vecz[] = { 0.025,-0.025,0.05,-0.025,0.05,-0.025,0.025,-0.025,0.025,-0.025,-0.025,-0.025,-0.025,-0.025,0.05,0.025,0.05,0.025 };

	glBegin(GL_QUAD_STRIP);
	for (int i = 0; i < num_of_vectors; i++) {
		setVertex(vecx[i], vecy[i], vecz[i]);
	}
	glEnd();

	//cannon gun
	drawGun(120, 0.005, 0.09);

	glRotatef(360 - rad2deg(dydz), 1, 0, 0);
	glRotatef(360 - rad2deg(dydx), 0, 0, 1);
	glRotatef(360-heading, 0, 1, 0);
	
	glTranslatef(-x, -y, -z);

	glColor3f(1.0, 1.0, 1.0);

}

void drawHealthBar()
{
	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	//glRotatef(45 - c.x, 0, 1, 0);
	//glRotatef(360 - c.y, 1, 0, 0);

	float displayHealthValue = health * 0.01 * 0.5;

	glTranslatef(-0.9 + (displayHealthValue / 2), 0.9, 0);
	glScalef(displayHealthValue, 0.04, 0.0);
	if (health > 40) glColor3f(0.0, 1.0, 0.0);
	else if (health > 20) glColor3f(1.0, 1.0, 0.0);
	else glColor3f(1.0 * (sinf((int)(glutGet(GLUT_ELAPSED_TIME) * 0.01)) * 0.5 + 0.5), 0.0, 0.0);
	//glColor3f(0.0, 1.0, 0.0);

	glBegin(GL_QUADS);
	glVertex2f(-0.5, -0.5);
	glVertex2f(0.5, -0.5);
	glVertex2f(0.5, 0.5);
	glVertex2f(-0.5, 0.5);
	glEnd();

	glColor3f(1.0, 1.0, 1.0);

	/* Pop modelview */
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);

	/* Pop projection */
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	/* Pop attributes */
	glPopAttrib();
}


void displayOSD()
{
	char buffer[30];
	char* bufp;
	int w, h;

	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	/* Set up orthographic coordinate system to match the
	   window, i.e. (0,0)-(w,h) */
	w = glutGet(GLUT_WINDOW_WIDTH);
	h = glutGet(GLUT_WINDOW_HEIGHT);
	glOrtho(0.0, w, 0.0, h, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	// Score points
	glColor3f(1.0, 1.0, 0.0);
	glRasterPos2i(20, h - 40);
	snprintf(buffer, sizeof buffer, "Score: %d", score);
	for (bufp = buffer; *bufp; bufp++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *bufp);

	/* Frame rate */
	glColor3f(1.0, 1.0, 0.0);
	glRasterPos2i(w-160, h-20);
	snprintf(buffer, sizeof buffer, "fr (f/s): %6.0f", global.frameRate);
	for (bufp = buffer; *bufp; bufp++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *bufp);

	/* Time per frame */
	glColor3f(1.0, 1.0, 0.0);
	glRasterPos2i(w - 160, h - 40);
	snprintf(buffer, sizeof buffer, "ft (ms/f): %5.0f", 1.0 / global.frameRate * 1000.0);
	for (bufp = buffer; *bufp; bufp++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *bufp);

	/* Tessalation */
	glColor3f(1.0, 1.0, 0.0);
	glRasterPos2i(w - 160, h - 60);
	snprintf(buffer, sizeof buffer, "tess: %i", n);
	for (bufp = buffer; *bufp; bufp++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *bufp);

	/* Pop modelview */
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);

	/* Pop projection */
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	/* Pop attributes */
	glPopAttrib();
}

void displayGameOver()
{
	char buffer[30];
	char* bufp;
	int w, h;

	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	/* Set up orthographic coordinate system to match the
	   window, i.e. (0,0)-(w,h) */
	w = glutGet(GLUT_WINDOW_WIDTH);
	h = glutGet(GLUT_WINDOW_HEIGHT);
	glOrtho(0.0, w, 0.0, h, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	/* Frame rate */
	glColor3f(1.0, 1.0, 1.0);
	glRasterPos2i(20, h-20);
	snprintf(buffer, sizeof buffer, "GAME OVER");
	for (bufp = buffer; *bufp; bufp++)
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *bufp);

	/* Pop modelview */
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);

	/* Pop projection */
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	/* Pop attributes */
	glPopAttrib();
}

void display(void)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	cameraFunc();

	drawSkyBox();

	if (wireframe) {
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		drawAxes(10.0f); //Main axes
		displayVelocity(cannon_angle, 0);
	}
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		if (lightview) {
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glEnable(GL_NORMALIZE);
		}
		else {
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
		}
	}

	drawHealthBar();

	glColor3f(1.0, 1.0, 1.0);

	drawGrid();

	drawIsland();

	drawTrajectoryAnalytical(cannon_angle, 0.0, 0.0);

	
	for (int i = 0; i < boatCount; i++)
	{
		if (!enemies[i].hit) drawBoat(enemies[i].x, enemies[i].z);
	}

	if (!gameOver) drawBall();

	displayOSD();

	if (gameOver) displayGameOver();

	glLoadIdentity();

	/* Always check for errors! */
	int err;
	while ((err = glGetError()) != GL_NO_ERROR)
		printf("display: %s\n", gluErrorString(err));

	glutSwapBuffers();

	//update framerate info
	global.frames++;
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'p':
		wireframe = !wireframe;
		printf("%d", wireframe);
		break;
	case 'n':
		axes = !axes;
		normalview = !normalview;
		tangentview = !tangentview;
		break;
	case 't':
		textureview = !textureview;
		break;
	case 'l':
		lightview = !lightview;
		break;
	
	case ' ':
		//fire cannon
		cannonballs[cnnbllCount] = {
		  { 0.0, 0.0, 0.0 },
		  { 0.25, 0.25, 0.25 },
		  { 0.0, 0.0, 0.0 },
		  {  0.25, 0.25, 0.0 },
		  cannon_angle.x,
		  cannon_angle.y,
		  waveTime.t,
		  initVel.speed,
		  FALSE
		};
		cnnbllCount++;
		break;


	case 'w':
		wPressed = TRUE;
		break;
	case 's':
		sPressed = TRUE;
		break;
	case 'a':
		aPressed = TRUE;
		break;
	case 'd':
		dPressed = TRUE;
		break;

	//tesselation control
	case '+':
	case '=':
		if (n < 1028) n *= 2;
		break;
	case '-':
		if (n > 2) n /= 2;
		break;
	case 27:
	case 'q':
		exit(EXIT_SUCCESS);
		break;
	default:
		break;
	}
}

void keyboardUp(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'w':
		wPressed = FALSE;
		break;
	case 'a':
		aPressed = FALSE;
		break;
	case 's':
		sPressed = FALSE;
		break;
	case 'd':
		dPressed = FALSE;
		break;
	}
}

void SpecialInput(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		break;
	case GLUT_KEY_DOWN:
		break;
	case GLUT_KEY_LEFT:
		leftKeyPressed = TRUE;
		break;
	case GLUT_KEY_RIGHT:
		rightKeyPressed = TRUE;
		break;
	}
}

void specialInputUp(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		break;
	case GLUT_KEY_DOWN:
		break;
	case GLUT_KEY_LEFT:
		leftKeyPressed = FALSE;
		break;
	case GLUT_KEY_RIGHT:
		rightKeyPressed = FALSE;
		break;
	}
}

void mouse(int button, int state, int x, int y)
{
	//store intial position when clicked
	if (button == GLUT_LEFT_BUTTON)
	{
		mouseMode = 0;
		c.x0 = x;
		c.y0 = y;
		//c.y = y;
	}

	if (button == GLUT_RIGHT_BUTTON)
	{
		mouseMode = 1;
		mouseZoom.x = x;
		mouseZoom.y = y;
	}

}


void mouseMotion(int x, int y)
{
	
	if (mouseMode == 0)
	{
		// update camera x angle based on last position
		c.x += x - c.x0;
		c.x0 = x;
		//limit camera y angle
		if ((c.y + y - c.y0) < 90 && (c.y + y - c.y0) > 5) c.y += y - c.y0;
		c.y0 = y;
	}
	else
	{
		//calculate and update camera zoom
		float dx = x - mouseZoom.x;
		mouseZoom.x = x;
		float dy = y - mouseZoom.y;
		mouseZoom.y = y;
		if ((c.zoom - (dx + dy) * 0.01) > -1 && (c.zoom - (dx + dy) * 0.01) < 0) {
			c.zoom -= (dx + dy) * 0.01;
		}
	}

}

boat spawnPosition()
{
	//generate a random boat spwan point
	/* initialize random seed: */
	srand(glutGet(GLUT_ELAPSED_TIME));
	float x = 0;
	float z = 0;

	int switcher = rand() % 2;

	/* generate decimal value of location*/
	if (switcher == 0) {
		while (x < 0.8 && x > -0.8) {
			x = rand() % 20;
			x = (float)x * 0.1 - 1;
		}
		z = rand() % 20;
		z = (float)z * 0.1 - 1;
	}
	else {
		while (z < 0.8 && z > -0.8) {
			z = rand() % 20;
			z = (float)z * 0.1 - 1;
		}
		x = rand() % 20;
		x = (float)x * 0.1 - 1;
	}
	
	
	return { x, z, FALSE };
}

void moveEnemies()
{
	for (int i = 0; i < boatCount; i++)
	{
		if (!enemies[i].hit) {
			//calculate distance between the boat and island
			float x = enemies[i].x;
			float z = enemies[i].z;
			float targetDistance = sqrt(pow(x, 2) + pow(z, 2));

			if (targetDistance >= 0.5) {
				if (enemies[i].x != 0) enemies[i].x < 0 ? enemies[i].x += 0.0005 : enemies[i].x -= 0.0005;
				if (enemies[i].z != 0) enemies[i].z < 0 ? enemies[i].z += 0.0005 : enemies[i].z -= 0.0005;
			}
			else {
				//calculate boat heading and make boats move in circle by increasing heading values;
				float heading = 0;
				heading = atanf(x/z) + 0.003;
				if (z < 0) heading -= M_PI;

				enemies[i].x = 0.5 * sinf(heading);
				enemies[i].z = 0.5 * cosf(heading);
			}
		}
	}
}

GLboolean detectCollision(state cnnbll)
{
	GLboolean collision = FALSE;
	float dydx = 0;
	float dydz = 0;

	//find height of sine wave, y
	float y = 0;
	float sum = 0;
	for (int wave = 0; wave < 2; wave++) {
		calcSineWave(sw[wave], cnnbll.r.x, cnnbll.r.z, waveTime.t, &y, TRUE, &dydx, &dydz);
		sum += y;
	}
	y = sum;

	if (cnnbll.r.y < y) collision = TRUE;

	return collision;

}

void detectHitBox(state cnnbll)
{
	float x = cnnbll.r.x;
	float z = cnnbll.r.z;
	//if the cnnbll travel shorter than island radius, island is hit
	if (sqrt(pow(x, 2) + pow(z, 2)) < island_r) {
		printf("Island hit\n");
		health--;
		printf("Health: %d\n", health);
	}
	else {
		//check which boat the cnnbll hit
		for (int i = 0; i < boatCount; i++) {
			if (!enemies[i].hit) {
				float ex = enemies[i].x;
				float ez = enemies[i].z;
				if ((x < ex + 0.05 && x > ex - 0.05) && (z < ez + 0.075 && z > ez - 0.075)) {
					printf("Boat hit\n");
					enemies[i].hit = TRUE;
					score += 100;
				}
			}
		}
	}

}

void calcCannonball()
{
	for (int i = 0; i < cnnbllCount; i++) {
		state cnnbll = cannonballs[i];
		//printf("cnnball phi");
		if (!cnnbll.sunk) {
			if (detectCollision(cnnbll) && !(cnnbll.r.x == 0 && cnnbll.r.z == 0)) {
				detectHitBox(cnnbll);
				cannonballs[i].sunk = TRUE;
			}
			else {
				updateProjectileStateAnalytical(waveTime.t - cnnbll.launchedTime, cnnbll, i, 0, cnnbll.cnnPower);
			}
		}
	}
}

void fireAtIsland()
{
	for (int i = 0; i < boatCount; i++) {
		
		if (!enemies[i].hit) {

			//generate a random chance to fire
			// initialize random seed:
			srand(glutGet(GLUT_ELAPSED_TIME) + i);
			int chance = rand() %3;

			if (chance == 1) {
				float ex = enemies[i].x;
				float ez = enemies[i].z;

				//calculate the angle required to fire at the island
				float d = sqrt(pow(ex, 2) + pow(ez, 2)) - 0.15;
				float power = boatVel.speed;
				printf("boat distance = %.2f", d);
				float angle = 90 - rad2deg(0.5 * asinf((-g * d) / pow(power, 2)));
				printf("boat angle = %.2f", angle);

				//calculate heading
				float heading = rad2deg(atan(ex / ez)) - 90;
				if (ez >= 0) heading -= 180;

				//find height of sine wave
				float ey = 0;
				float sum = 0;
				float dydx = 0;
				float dydz = 0;

				for (int wave = 0; wave < 2; wave++) {
					calcSineWave(sw[wave], ex, ez, waveTime.t, &ey, FALSE, &dydx, &dydz);
					sum += ey;
				}
				ey = sum + 0.001;

				cannonballs[cnnbllCount] = {
				  { ex, ey, ez },
				  { 0.25, 0.25, 0.25 },
				  { ex, ey, ez },
				  {  0.25, 0.25, 0.0 },
				  angle,
				  heading,
				  waveTime.t,
				  power,
				  FALSE
				};
				cnnbllCount++;
				printf("boat fired their cannon\n");
			}
		}
	}
}

void processKeyboard()
{
	if (wPressed) initVel.speed += 0.0015;
	if (sPressed) initVel.speed -= 0.0015;
	if (aPressed && cannon_angle.x <= 90) cannon_angle.x+=0.1;
	if (dPressed && cannon_angle.x >= 0) cannon_angle.x-=0.1;
	if (leftKeyPressed) { 
		cannon_angle.y += 0.5;
		c.x-=0.5; 
	}
	if (rightKeyPressed) {
		cannon_angle.y -= 0.5;
		c.x+=0.5;
	}
}

void idle()
{
	float t, dt = 0;

	/* Frame rate */
	t = glutGet(GLUT_ELAPSED_TIME) / (float)milli - global.startTime;

	if (lastT < 0.0) {
		lastT = t;
		return;
	}

	dt = t - lastT;

	if (!gameOver) {
		processKeyboard();

		//animate waves
		waveTime.t = glutGet(GLUT_ELAPSED_TIME) * 0.001;

		//Create and spawn enemy boats
		if (waveTime.t > boatSpawnTimer&& boatCount < (int)(sizeof(enemies) / sizeof(enemies[0])))
		{
			enemies[boatCount] = spawnPosition();
			boatCount++;
			boatSpawnTimer += 3;	//3 second timer
		}

		//make enemy boats fire at island
		if (waveTime.t > shootTimer) {
			fireAtIsland();
			glutPostRedisplay();
			shootTimer += 1;
		}

		//cannonball
		calcCannonball();

		glutPostRedisplay();

		moveEnemies();

		if (health <= 0) gameOver = TRUE;

		glutPostRedisplay();
	}

	lastT = t;
	dt = t - global.lastFrameRateT;

	if (dt > global.frameRateInterval) {
		global.frameRate = global.frames / dt;
		global.lastFrameRateT = t;
		global.frames = 0;
	}

}

void init()
{
	/* In this program these OpenGL calls only need to be done once,
	  but normally they would go elsewhere, e.g. display */

	glMatrixMode(GL_PROJECTION);
	glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);

	global.startTime = glutGet(GLUT_ELAPSED_TIME) / (float)milli;
	//printf("global start time: %.2f", global.startTime);
	global.go = TRUE;
}

void reshape(int width, int height)
{
	glViewport(0, 0, width, height);

	aspect = ((float)width / (float)height);

	glMatrixMode(GL_PROJECTION);	// Set the matrix we want to manipulate
	glLoadIdentity();			// Overwrite the contents with the identity
	gluPerspective(fovy, aspect, zNear, zFar);		// Multiply the current matrix with a generated perspective matrix
	glMatrixMode(GL_MODELVIEW);	// Change back to the modelview matrix

	glutPostRedisplay();
}

int loadSkyBoxTexture()
{
	skyboxTexture.xn = loadTexture("textures/xneg.png");
	skyboxTexture.xp = loadTexture("textures/xpos.png");
	skyboxTexture.yn = loadTexture("textures/yneg.png");
	skyboxTexture.yp = loadTexture("textures/ypos.png");
	skyboxTexture.zn = loadTexture("textures/zneg.png");
	skyboxTexture.zp = loadTexture("textures/zpos.png");

	if (!skyboxTexture.xn || !skyboxTexture.xp || !skyboxTexture.yn || !skyboxTexture.yp || !skyboxTexture.zn || !skyboxTexture.zp) {
		printf("No texture created; exiting.\n");
		return 1;
	}

	return 0;
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("Pearl Harbour 2: Nuclear War");	//A funnier title compared to "I3D Assignment 2"

	init();

	if (!gameOver) {
		glutReshapeFunc(reshape);
		glutDisplayFunc(display);
		glutKeyboardFunc(keyboard);
		glutKeyboardUpFunc(keyboardUp);
		glutSpecialFunc(SpecialInput);
		glutSpecialUpFunc(specialInputUp);
		glutMotionFunc(mouseMotion);
		glutMouseFunc(mouse);
		glutIdleFunc(idle);

		texture = loadTexture("textures/sand.jpg");
		if (!texture) {
			printf("No texture created; exiting.\n");
			return EXIT_FAILURE;
		}

		int skybox = loadSkyBoxTexture();
		if (skybox == EXIT_FAILURE) {
			printf("No texture created; exiting.\n");
			return EXIT_FAILURE;
		}

		glutMainLoop();
	}

	return EXIT_SUCCESS;
}
