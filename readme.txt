COSC1187 Interactive 3D Graphics and Animation
Assignment 2: Island Defence 3D (a.k.a Pearl Harbour 2: Nuclear War)

Name: Raphael Doh Ong Wong
ID: s3735236

-------------
Description
-------------
This project is built for assignment 2 for COSC1187 Interactive 3D Graphics and Animation. This program is built for education
purposes only.

Assignment 2: Island Defence 3D (a.k.a Pearl Harbour 2: Nuclear War)
--------------------------------------------------------------------

This is an island survival game where enemy battleships attack and surrounds you. Defend yourself with a powerful cannon and sunk your 
enemies! Show no mercy!

-------------------
Installation guide:
-------------------
Setup Xming server (Windows)
- Run Xming
- type 'export DISPLAY=:0.0' in WSL

Compile:
"make"
-or-
"g++ -Wall -o a2 *.cpp -lglut -lGLU -lGL -lm"

Run:
./a2

-----------
Features
-----------
- Animated Waves now in 3D!
- Advanced AI enemy
- OSD display
- Cannon ball shoots from the end of the cannon (main island only)
- Wave tesselation control
- Skybox
- HD Sand Texture
- Colourful Health Bar

--------
Controls
--------

Player:
	left 		rotate cannon to left
	right		rotate cannon to right
	mouse		camera control
	left click	enable camera control
	right click	camera zoom
	w		increase cannon power
	s		decrease cannon power
	a		increase cannon angle
	d		decrease cannon angle
	space		fire cannon

Debug:
	p		wireframe mode
	n		normal and tanget view
	t		texture
	l		light

------------
Known / Potential Issues
------------
- Opaque waters
- Drops in frame rate when window is resized due to resolution increase
- Speed of control changes depending on resolution / window size
- MVC not implemented (this programme was developped on Visual Studio and I cannot figure how to use project explorer.)
- Shininess not working
- Makefile is tested in Ubuntu for Windows
- No restart option after game over
- Normal / tangent vision only applies for waves

