
#include <string.h>
#include "camera.h"
//#include "Shape.h"

unsigned char frame[windowW * windowH * 3];

void OnKeyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'q':
		exit(0);
		break;
	}
}

void DrawPixel(unsigned int x, unsigned int y, unsigned char r, unsigned char g, unsigned char b)
{
	if (x >= windowW || y >= windowH)
		return;

	unsigned int index = 3 * (y * windowW + x);
	frame[index] = r;
	frame[index + 1] = g;
	frame[index + 2] = b;
}

void Bresenham(int x1, int y1, int x2, int y2, unsigned char r, unsigned char g, unsigned char b)
{

	// To remark end points as red.
	//    DrawPixel(x1,y1,255,0,0);
	//    DrawPixel(x1+1,y1,255,0,0);
	//    DrawPixel(x1,y1+1,255,0,0);
	//    DrawPixel(x1+1,y1+1,255,0,0);
	//
	//    DrawPixel(x2,y2,255,0,0);
	//    DrawPixel(x2+1,y2,255,0,0);
	//    DrawPixel(x2,y2+1,255,0,0);
	//    DrawPixel(x2+1,y2+1,255,0,0);

	int dx = x2 - x1,
		dy = y2 - y1,
		y = y1,
		epsilon = 0; // epsilon as error.

	/*
	 * As we know, a Bresenham's Algo must work together for 8 octants. Starts from 12 o' clock, I named the first
	 * octants as #1, then #2; the 3'o clock's as #3, #4, lastly end with #8.
	 *
	 * Obviously, Case 2 and Case 1 is the basis of all cases. Case 2 could generate #3, #7, #6 and Case 1 do the rest.
	 *
	 * Case 2, I learn from class note, and I do the rest all by myself.
	 */

	 /* First four groups: based on Case 2 */

	if (0 <= dy && dy <= dx) { // 0 < m < 1, case 2
		//std::cout << "case 2" << std::endl;
		for (int x = x1; x <= x2; x += 1) {
			DrawPixel(x, y, r, g, b);
			epsilon += dy;
			if ((epsilon << 1) >= dx) { //Some little optimization, avoiding division with left-shift.
				y++;  epsilon -= dx;
			}
		}
	}

	if (dy >= 0 && dy <= -dx) { // 0 > m > -1, case 7

		// Case 7 is mirror of 2, around y-axis.  So change For Loop's order, and dx, dy into abs.

		//std::cout << "case 7" << std::endl;
		dx = -dx;
		for (int x = x1; x >= x2; x -= 1) {
			DrawPixel(x, y, r, g, b);
			epsilon += dy;
			if ((epsilon << 1) >= dx) {
				y++;  epsilon -= dx;
			}
		}
	}

	else if (dy <= 0 && -dy <= dx) { // -1 < m < 0, case 3

		// 3 is mirror of 2, around x-axis. So y++ => y--, and dx, dy into abs.

		//std::cout << "case 3" << std::endl;
		dy = -dy;
		for (int x = x1; x <= x2; x += 1) {
			DrawPixel(x, y, r, g, b);
			epsilon += dy;
			if ((epsilon << 1) >= dx) {
				y--;  epsilon -= dx;
			}
		}
	}

	//Another version of case 3 is regard it as case 7's reverse. So exchange the end points.
	//Bresenham(x2,y2, x1, y1);

	else if (dy <= 0 && dx <= dy) { // -1 < m < 0, case 6
		//case 6 is just a reverse of case 2. So exchange the end points.
		Bresenham(x2, y2, x1, y1, r, g, b);
	}


	/* Another four groups: based on Case 1
	 *
	 * The design idea of Case 1 is, it's a mirror around line y = x.
	 *
	 * So simply exchange x,y, and all is done.
	 * However, just call Bre(y1,x1,y2,x2) won't work for it.
	 */

	else if (0 <= dx && dx <= dy) { // m > 1, case 1

		//Case 1 is case 2's sibling: Change every x and y, but don't change Draw(x,y).

		//std::cout << "case 1" << std::endl;

		int x = x1;
		for (int y = y1; y <= y2; y += 1) {
			DrawPixel(x, y, r, g, b);
			epsilon += dx;
			if ((epsilon << 1) >= dy) {
				x++;  epsilon -= dy;
			}
		}
	}

	else if (dx <= 0 && -dx <= dy) { // m < -1, case 8

		//Case 8 is case 1's mirror around y-axis. Simply change the order of for loop; Abs dx, dy; and x = x1 => x = x2;

		//std::cout << "case 8" << std::endl;

		dx = -dx;
		int x = x2;
		for (int y = y2; y >= y1; y -= 1) {
			DrawPixel(x, y, r, g, b);
			epsilon += dx;
			if ((epsilon << 1) >= dy) {
				x++;  epsilon -= dy;
			}
		}
	}

	else if (dx <= 0 && dy <= dx) { // m < -1, case 5
		//case 5 is just a reverse of case 1. So exchange the end points.
		Bresenham(x2, y2, x1, y1, r, g, b);
	}

	else if (0 <= dx && dy <= -dx) { // m < -1, case 4
		//case 4 is just a reverse of case 8. So exchange the end points.
		Bresenham(x2, y2, x1, y1, r, g, b);
	}
}

/*
   Draw a line by Bresenham's algorithm.

   @author Mu He

   This function is writen because an error in our assignment.
   In assignment#1 the bresenham algo ask us input int, but know is double/float

   @author Mu He
  */
void Line(double x1, double y1, double x2, double y2)
{
	Bresenham((int)x1, (int)y1, (int)x2, (int)y2, 0, 0, 0);
}

void torus()
{
	dmatrix_t x, y, z;//x, y, z for torus

	dmat_alloc(&x, 4, 1);
	dmat_alloc(&y, 4, 1);
	dmat_alloc(&z, 4, 1);
	//parametric equation
		//loc
	double radium_torus, radium_torus_pipe;
	radium_torus = 200; radium_torus_pipe = 20;

	double du = 2.0*M_PI / 8;
	double dv = 2.0*M_PI / 100;

	for (double v = 0.0; v <= 2.0*M_PI; v += dv) {
		for (double u = 0.0; u <= 2.0*M_PI; u += du) {
			//x, y, z for torus
			x.m[1][1] = (radium_torus + radium_torus_pipe * cos(u))*cos(v);
			x.m[2][1] = (radium_torus + radium_torus_pipe * cos(u))*sin(v);
			x.m[3][1] = radium_torus_pipe * sin(u);
			x.m[4][1] = 1;
			x = *perspective_projection(dmat_mult(&C, &x));
			y.m[1][1] = (radium_torus + radium_torus_pipe * cos(u + du))*cos(v);
			y.m[2][1] = (radium_torus + radium_torus_pipe * cos(u + du))*sin(v);
			y.m[3][1] = radium_torus_pipe * sin(u + du);
			y.m[4][1] = 1;
			y = *perspective_projection(dmat_mult(&C, &y));
			z.m[1][1] = (radium_torus + radium_torus_pipe * cos(u))*cos(v + dv);
			z.m[2][1] = (radium_torus + radium_torus_pipe * cos(u))*sin(v + dv);
			z.m[3][1] = radium_torus_pipe * sin(u);
			z.m[4][1] = 1;
			z = *perspective_projection(dmat_mult(&C, &z));

			Line(x.m[1][1], x.m[2][1], y.m[1][1], y.m[2][1]);
			Line(x.m[1][1], x.m[2][1], z.m[1][1], z.m[2][1]);
		}
	}

	delete_dmatrix(&x);
	delete_dmatrix(&y);
	delete_dmatrix(&z);
}

void sphere()
{
	double radium_shpere = 80;

	dmatrix_t x_s, y_s, z_s;//x_s, y_s, z_s for shpere
	dmat_alloc(&x_s, 4, 1);
	dmat_alloc(&y_s, 4, 1);
	dmat_alloc(&z_s, 4, 1);

	//u_s, v_s, du_s, dv_s are values for shpere
	double dv_s = 2.0*M_PI / 40.0;
	double du_s = M_PI / 16.0;

	/* Notice: the shpere's u, that is every slices' euqation, only need run a half */
	for (double v_s = 0.0; v_s <= 2 * M_PI; v_s += dv_s) {
		for (double u_s = 0.0; u_s <= M_PI; u_s += du_s) {
			//loc
			x_s.m[1][1] = radium_shpere * cos(v_s)*sin(u_s);
			x_s.m[2][1] = radium_shpere * sin(v_s)*sin(u_s);
			x_s.m[3][1] = radium_shpere * cos(u_s);
			x_s.m[4][1] = 1;
			x_s = *perspective_projection(dmat_mult(&C, &x_s));
			y_s.m[1][1] = radium_shpere * cos(v_s)*sin(u_s + du_s);
			y_s.m[2][1] = radium_shpere * sin(v_s)*sin(u_s + du_s);
			y_s.m[3][1] = radium_shpere * cos(u_s + du_s);
			y_s.m[4][1] = 1;
			y_s = *perspective_projection(dmat_mult(&C, &y_s));
			z_s.m[1][1] = radium_shpere * cos(v_s + dv_s)*sin(u_s);
			z_s.m[2][1] = radium_shpere * sin(v_s + dv_s)*sin(u_s);
			z_s.m[3][1] = radium_shpere * cos(u_s);
			z_s.m[4][1] = 1;
			z_s = *perspective_projection(dmat_mult(&C, &z_s));

			Line(x_s.m[1][1],
				x_s.m[2][1],
				y_s.m[1][1],
				y_s.m[2][1]);
			Line(x_s.m[1][1],
				x_s.m[2][1],
				z_s.m[1][1],
				z_s.m[2][1]);
		}
	}
	delete_dmatrix(&x_s);
	delete_dmatrix(&y_s);
	delete_dmatrix(&z_s);
}

void Draw()
{
    torus(); sphere();
//    Torus t;
//    Sphere s;
//
//    t.draw(1, 1, 1, 1);
//    s.draw(1, 1, 1);
}

void OnDisplay()
{
	memset(frame, 255, windowW * windowH * 3);
	//main print function
	Draw();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(windowW, windowH, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte*)frame);
    glFlush();
    glutSwapBuffers();
}
