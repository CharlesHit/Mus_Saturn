#include <vector>
#include <X11/Xlib.h>
#include <OpenGL/gl.h>
#include <OpenGl/glu.h>
#include <GLUT/glut.h>
#include <string.h>
#include "camera.h"
using namespace std;

int REMOVAL = 1;
const int X = 1;
const int Y = 2;

double angle_total = 0.0;  //totally rotated angle

Display *d;
Window w;
XEvent e;
int s;

vector<char> frame(windowW * windowH * 3);
vector<double> depth(windowW * windowH);

dmatrix_t *crossProduct(dmatrix_t vect_A, dmatrix_t vect_B)
{
    dmatrix_t result;
    dmat_alloc(&result, 4, 1);
    result.m[1][1] = vect_A.m[2][1] * vect_B.m[3][1] - vect_A.m[3][1] * vect_B.m[2][1];
    result.m[2][1] = vect_A.m[1][1] * vect_B.m[3][1] - vect_A.m[3][1] * vect_B.m[1][1];
    result.m[3][1] = vect_A.m[1][1] * vect_B.m[2][1] - vect_A.m[2][1] * vect_B.m[1][1];
    result.m[4][1] = 0;

    return &result;
}

double angle(dmatrix_t vect_A, dmatrix_t vect_B)
{
    dmatrix_t vecLight_xyz;
    dmat_alloc(&vecLight_xyz, 3, 1);
    vecLight_xyz.m[1][1] = vect_A.m[1][1];
    vecLight_xyz.m[2][1] = vect_A.m[2][1];
    vecLight_xyz.m[3][1] = vect_A.m[3][1];

    dmatrix_t vecNorm_xyz;
    dmat_alloc(&vecNorm_xyz, 3, 1);
    vecNorm_xyz.m[1][1] = vect_B.m[1][1];
    vecNorm_xyz.m[2][1] = vect_B.m[2][1];
    vecNorm_xyz.m[3][1] = vect_B.m[3][1];

    return ddot_product(&vecLight_xyz, &vecNorm_xyz) / ((dmat_norm(&vecLight_xyz)) * (dmat_norm(&vecNorm_xyz)));
}

void OnKeyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 'q':
            printf("total rotate angle: %f\n", angle_total);
            exit(0);
            break;
    }
}

void DrawPixel(unsigned int x, unsigned int y, unsigned char r, unsigned char g, unsigned char b, double floatDepth)
{
    if (x >= windowW || y >= windowH)
        return;

    long long index = (y * windowW + x);
    if (floatDepth < depth[index])
    {
        depth[index] = floatDepth;

        index = 3 * index;
        frame[index] = r;
        frame[index + 1] = g;
        frame[index + 2] = b;
    }
}

void Line(int x1, int y1, int x2, int y2, unsigned char r, unsigned char g, unsigned char b, double floatDepth_start,
          double floatDepth_end)
{
    double floatDepth = floatDepth_end - floatDepth_start;
    double dDepth;
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

    if (0 <= dy && dy <= dx)
    { // 0 < m < 1, case 2
        //std::cout << "case 2" << std::endl;

        dDepth = floatDepth / (x2 - x1);

        for (int x = x1; x <= x2; x += 1)
        {
            DrawPixel(x, y, r, g, b, floatDepth_start + dDepth);
            epsilon += dy;
            if ((epsilon << 1) >= dx)
            { //Some little optimization, avoiding division with left-shift.
                y++;
                epsilon -= dx;
            }
        }
    }

    if (dy >= 0 && dy <= -dx)
    { // 0 > m > -1, case 7

        // Case 7 is mirror of 2, around y-axis.  So change For Loop's order, and dx, dy into abs.

        //std::cout << "case 7" << std::endl;
        dx = -dx;

        dDepth = floatDepth / (x2 - x1);

        for (int x = x1; x >= x2; x -= 1)
        {
            DrawPixel(x, y, r, g, b, floatDepth_start + dDepth);
            epsilon += dy;
            if ((epsilon << 1) >= dx)
            {
                y++;
                epsilon -= dx;
            }
        }
    } else if (dy <= 0 && -dy <= dx)
    { // -1 < m < 0, case 3

        // 3 is mirror of 2, around x-axis. So y++ => y--, and dx, dy into abs.

        //std::cout << "case 3" << std::endl;
        dy = -dy;

        dDepth = floatDepth / (x2 - x1);

        for (int x = x1; x <= x2; x += 1)
        {
            DrawPixel(x, y, r, g, b, floatDepth_start + dDepth);
            epsilon += dy;
            if ((epsilon << 1) >= dx)
            {
                y--;
                epsilon -= dx;
            }
        }
    }

        //Another version of case 3 is regard it as case 7's reverse. So exchange the end points.
        //Bresenham(x2,y2, x1, y1);

    else if (dy <= 0 && dx <= dy)
    { // -1 < m < 0, case 6
        //case 6 is just a reverse of case 2. So exchange the end points.
        Line(x2, y2, x1, y1, r, g, b, floatDepth_end, floatDepth_start);
    }


        /* Another four groups: based on Case 1
         *
         * The design idea of Case 1 is, it's a mirror around line y = x.
         *
         * So simply exchange x,y, and all is done.
         * However, just call Bre(y1,x1,y2,x2) won't work for it.
         */

    else if (0 <= dx && dx <= dy)
    { // m > 1, case 1

        //Case 1 is case 2's sibling: Change every x and y, but don't change Draw(x,y).

        //std::cout << "case 1" << std::endl;

        dDepth = floatDepth / (x2 - x1);

        int x = x1;
        for (int y = y1; y <= y2; y += 1)
        {
            DrawPixel(x, y, r, g, b, floatDepth_start + dDepth);
            epsilon += dx;
            if ((epsilon << 1) >= dy)
            {
                x++;
                epsilon -= dy;
            }
        }
    } else if (dx <= 0 && -dx <= dy)
    { // m < -1, case 8

        //Case 8 is case 1's mirror around y-axis. Simply change the order of for loop; Abs dx, dy; and x = x1 => x = x2;

        //std::cout << "case 8" << std::endl;

        dDepth = floatDepth / (x2 - x1);

        dx = -dx;
        int x = x2;
        for (int y = y2; y >= y1; y -= 1)
        {
            DrawPixel(x, y, r, g, b, floatDepth_start + dDepth);
            epsilon += dx;
            if ((epsilon << 1) >= dy)
            {
                x++;
                epsilon -= dy;
            }
        }
    } else if (dx <= 0 && dy <= dx)
    { // m < -1, case 5
        //case 5 is just a reverse of case 1. So exchange the end points.
        Line(x2, y2, x1, y1, r, g, b, floatDepth_end, floatDepth_start);
    } else if (0 <= dx && dy <= -dx)
    { // m < -1, case 4
        //case 4 is just a reverse of case 8. So exchange the end points.
        Line(x2, y2, x1, y1, r, g, b, floatDepth_end, floatDepth_start);
    }
}

/*
   Draw a line by Bresenham's algorithm.

   @author Mu He

   This function is writen because an error in our assignment.
   In assignment#1 the bresenham algo ask us input int, but know is double/float

   @author Mu He
  */
void Line(double x1, double y1, double x2, double y2, double floatDepth_start, double floatDepth_end)
{
    Line((int) x1, (int) y1, (int) x2, (int) y2, 0, 0, 0, floatDepth_start, floatDepth_end);
}

void Line3D(dmatrix_t x, dmatrix_t y, dmatrix_t matProj)
{
    double floatDepth_start, floatDepth_end;
    floatDepth_start = dmat_norm(dmat_sub(&x, &pointEye));
    floatDepth_end = dmat_norm(dmat_sub(&y, &pointEye));

    x = *perspective_projection(dmat_mult(&matProj, &x));
    y = *perspective_projection(dmat_mult(&matProj, &y));

    Line((int) x.m[1][1], (int) x.m[2][1], (int) y.m[1][1], (int) y.m[2][1], 0, 0, 0, floatDepth_start, floatDepth_end);
}

double minimum_coordinate(int coordinate, dmatrix_t P[], int n)
{

    int i;
    double min;

    min = P[0].m[coordinate][1];

    for (i = 1; i < n; i++)
    {
        if (P[i].m[coordinate][1] < min)
        {
            min = P[i].m[coordinate][1];
        }
    }
    return min;
}

double maximum_coordinate(int coordinate, dmatrix_t P[], int n)
{

    int i;
    double max;

    max = P[0].m[coordinate][1];

    for (i = 1; i < n; i++)
    {
        if (P[i].m[coordinate][1] > max)
        {
            max = P[i].m[coordinate][1];
        }
    }
    return max;
}

int maximum_intersection(int intersections[], int n)
{

    int i, max;

    max = intersections[0];
    for (i = 1; i < n; i++)
    {
        if (intersections[i] > max)
        {
            max = intersections[i];
        }
    }
    return max;
}

int minimum_intersection(int intersections[], int n)
{

    int i, min;

    min = intersections[0];
    for (i = 1; i < n; i++)
    {
        if (intersections[i] < min)
        {
            min = intersections[i];
        }
    }
    return min;
}

void
XFillConvexPolygon(Display *d, Window w, int s, dmatrix_t xx, dmatrix_t yy, dmatrix_t zz, int n, unsigned int r, unsigned int g, unsigned int bb,
                   double floatDepth)
{
    dmatrix_t P[3];
    for (int i = 0; i < 3; i++) dmat_alloc(&P[i], 4, 1);
    P[0] = xx;
    P[1] = yy;
    P[2] = zz;

    unsigned int i, j;
    int y, y_min, y_max, min_int, max_int;
    double m, b;
    int *active, *horizontal, *intersections;

    horizontal = (int *) malloc(n * sizeof(int)); /* Allocate horizontal segment table */
    active = (int *) malloc(n * sizeof(int)); /* Allocate active segment table */
    intersections = (int *) malloc(n * sizeof(int)); /* Allocate intersection table */

    y_min = (int) minimum_coordinate(Y, P, n); /* Determine number of scan lines */
    y_max = (int) maximum_coordinate(Y, P, n);

    for (i = 0; i < n; i++)
    {
        horizontal[i] = (int) P[i].m[Y][1] == (int) P[(i + 1) % n].m[Y][1]; /* Find horizontal segments */
    }

    for (y = y_min; y <= y_max; y++)
    { /* For each scan line y */
        for (i = 0; i < n; i++)
        {  /* Update segment table */
            if (!horizontal[i])
            {
                active[i] = (y >= (int) P[i].m[Y][1] && y <= (int) P[(i + 1) % n].m[Y][1]) ||
                            (y <= (int) P[i].m[Y][1] && y >= (int) P[(i + 1) % n].m[Y][1]);
            }
        }
        j = 0;
        for (i = 0; i < n; i++)
        { /* find intersection x-value. The y-value is given by the scan line */
            if (active[i] && !horizontal[i])
            {
                if ((int) P[i].m[X][1] == (int) P[(i + 1) % n].m[X][1])
                { /* Vertical segment */
                    intersections[j++] = (int) P[i].m[X][1];
                } else
                {
                    m = (double) ((int) P[(i + 1) % n].m[Y][1] - (int) P[i].m[Y][1]) /
                        (double) ((int) P[(i + 1) % n].m[X][1] - (int) P[i].m[X][1]); /* Compute slope and intercept */
                    b = (double) ((int) P[i].m[Y][1]) - m * (double) ((int) P[i].m[X][1]);
                    intersections[j++] = (int) (((double) y - b) / m); /* Compute intersection */
                }
            }
        }
        min_int = minimum_intersection(intersections, j);
        max_int = maximum_intersection(intersections, j) + 1;
        for (i = min_int; i < max_int; i++)
        { /* Tracing from minimum to maximum intersection */
            DrawPixel(i, y, r, g, bb, floatDepth);
        }
    }
    free(horizontal);
    free(active);
    free(intersections);
}


//Notice: I use not classical way to get normal vector
void torus(unsigned int r, unsigned int g, unsigned int b)
{
    //information of torus
    double radium_torus, radium_torus_pipe;
    radium_torus = 180;
    radium_torus_pipe = 4;
    //the range of torus
    double v = 0.0 * M_PI;
    double u = 0.0 * M_PI;
    double du = 2.0 * M_PI / 8;
    double dv = 2.0 * M_PI / 64;
    double floatDepth_xyz, floatDepth_yzb, brightness_xyz, brightness_yzb;

    //initialization
    dmatrix_t x, y, z, buffer, center, vecYX, vecXZ, vecYZ, vecZbuffer, pointCenter_xyz, pointCenter_yzb, vecNorm_xyz, vecNorm_yzb, vecLight_xyz, vecLight_yzb, pointNormEnd_xyz, pointNormEnd_yzb;//x, y, z for torus, center of every circle on torus
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&vecYX, 4, 1);
    dmat_alloc(&vecXZ, 4, 1);
    dmat_alloc(&vecYZ, 4, 1);
    dmat_alloc(&vecZbuffer, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    dmat_alloc(&center, 4, 1);
    dmat_alloc(&pointCenter_xyz, 4, 1);
    dmat_alloc(&pointCenter_yzb, 4, 1);
    dmat_alloc(&vecNorm_xyz, 4, 1);
    dmat_alloc(&vecNorm_yzb, 4, 1);
    dmat_alloc(&vecLight_xyz, 4, 1);
    dmat_alloc(&vecLight_yzb, 4, 1);
    dmat_alloc(&pointNormEnd_xyz, 4, 1);
    dmat_alloc(&pointNormEnd_yzb, 4, 1);

    //the buffer is the x's opposite, so x, y, z, buffer makes a circle:
    // x - z
    // |
    // y - buffer
    buffer.m[1][1] = (radium_torus + radium_torus_pipe * cos(u)) * cos(v);
    buffer.m[2][1] = (radium_torus + radium_torus_pipe * cos(u)) * sin(v);
    buffer.m[3][1] = radium_torus_pipe * sin(u);;
    buffer.m[4][1] = 1;

    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = -100.0;
    pointLight.m[2][1] = -100.0;
    pointLight.m[3][1] = 300.0;
    pointLight.m[4][1] = 1.0;

    for (v = 0 * M_PI; v <= 2 * M_PI; v += dv)
    {
        for (u = 0 * M_PI; u <= 2 * M_PI; u += du)
        {
            //0.x, y, z calculation part;
            center.m[1][1] = radium_torus * cos(v);
            center.m[2][1] = radium_torus * sin(v);
            center.m[3][1] = 0;
            center.m[4][1] = 1;
            x.m[1][1] = (radium_torus + radium_torus_pipe * cos(u)) * cos(v);
            x.m[2][1] = (radium_torus + radium_torus_pipe * cos(u)) * sin(v);
            x.m[3][1] = radium_torus_pipe * sin(u);
            x.m[4][1] = 1;
            y.m[1][1] = (radium_torus + radium_torus_pipe * cos(u + du)) * cos(v);
            y.m[2][1] = (radium_torus + radium_torus_pipe * cos(u + du)) * sin(v);
            y.m[3][1] = radium_torus_pipe * sin(u + du);
            y.m[4][1] = 1;
            z.m[1][1] = (radium_torus + radium_torus_pipe * cos(u)) * cos(v + dv);
            z.m[2][1] = (radium_torus + radium_torus_pipe * cos(u)) * sin(v + dv);
            z.m[3][1] = radium_torus_pipe * sin(u);
            z.m[4][1] = 1;
            buffer.m[1][1] = (radium_torus + radium_torus_pipe * cos(u + du)) * cos(v + dv);
            buffer.m[2][1] = (radium_torus + radium_torus_pipe * cos(u + du)) * sin(v + dv);
            buffer.m[3][1] = radium_torus_pipe * sin(u + du);
            buffer.m[4][1] = 1;


            //1. edge vectors
            dmat_alloc(&vecYX, 4, 1);
            vecYX = *dmat_sub(&x, &y);//x-y
            vecXZ = *dmat_sub(&z, &x);
            vecYZ = *dmat_sub(&z, &y);
            vecZbuffer = *dmat_sub(&buffer, &z);


            //2.point center
            pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1]) / 3;
            pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1]) / 3;
            pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1]) / 3;
            pointCenter_xyz.m[4][1] = 1;

            pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1]) / 3;
            pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1]) / 3;
            pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1]) / 3;
            pointCenter_yzb.m[4][1] = 1;


            //3. Depth of each mesh
            floatDepth_xyz = dmat_norm(dmat_sub(&pointEye, &pointCenter_xyz));
            floatDepth_yzb = dmat_norm(dmat_sub(&pointEye, &pointCenter_yzb));


            //4.calculate Normal vectors
            //vecNorm_xyz = *crossProduct(vecYX, vecXZ); //the classical way
            vecNorm_xyz = *dmat_sub(&pointCenter_xyz, &center);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);
            vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //vecNorm_yzb = *crossProduct(vecYZ, vecZbuffer); //the classical way
            vecNorm_yzb = *dmat_sub(&pointCenter_yzb, &center);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);
            vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.


            //5. light vector
            vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);


            //6. brightness: cos of two vector X light
            brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            brightness_yzb = brightness_xyz; //angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;


            //7.light-center & normal vector (optional)
            //7.1 xyz
            pointNormEnd_xyz = *dmat_add(&pointCenter_xyz,
                                         &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1] ,0 ,0);
            //normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1],0,0);

            //7.2 yzb
            pointNormEnd_yzb = *dmat_add(&pointCenter_yzb,
                                         &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1] ,0 ,0);
            //normal vector
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1],0,0);


            //8. draw the meshes.
            //Line3D(x, y, C);
            //Line3D(x, z, C);
            //Line3D(y, z, C);


            //9.fill with color
            x = *perspective_projection(dmat_mult(&C, &x));
            y = *perspective_projection(dmat_mult(&C, &y));
            z = *perspective_projection(dmat_mult(&C, &z));
            buffer = *perspective_projection(dmat_mult(&C, &buffer));

            XFillConvexPolygon(d, w, s, x, y, z, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz, floatDepth_xyz);
            XFillConvexPolygon(d, w, s, buffer, y, z, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb, floatDepth_yzb);
        }
    }


    //9.clean
    free_dmatrix(x.m, 1, x.l, 1, x.c);
    free_dmatrix(y.m, 1, y.l, 1, y.c);
    free_dmatrix(z.m, 1, z.l, 1, z.c);
    free_dmatrix(buffer.m, 1, buffer.l, 1, buffer.c);
    free_dmatrix(pointLight.m, 1, pointLight.l, 1, pointLight.c);
    free_dmatrix(vecYX.m, 1, vecYX.l, 1, vecYX.c);
    free_dmatrix(vecXZ.m, 1, vecXZ.l, 1, vecXZ.c);
    free_dmatrix(vecYZ.m, 1, vecYZ.l, 1, vecYZ.c);
    free_dmatrix(vecZbuffer.m, 1, vecZbuffer.l, 1, vecZbuffer.c);
    free_dmatrix(center.m, 1, center.l, 1, center.c);
    free_dmatrix(pointCenter_xyz.m, 1, pointCenter_xyz.l, 1, pointCenter_xyz.c);
    free_dmatrix(pointCenter_yzb.m, 1, pointCenter_yzb.l, 1, pointCenter_yzb.c);
    free_dmatrix(vecNorm_xyz.m, 1, vecNorm_xyz.l, 1, vecNorm_xyz.c);
    free_dmatrix(vecNorm_yzb.m, 1, vecNorm_yzb.l, 1, vecNorm_yzb.c);
    free_dmatrix(vecLight_xyz.m, 1, vecLight_xyz.l, 1, vecLight_xyz.c);
    free_dmatrix(vecLight_yzb.m, 1, vecLight_yzb.l, 1, vecLight_yzb.c);
    free_dmatrix(pointNormEnd_xyz.m, 1, pointNormEnd_xyz.l, 1, pointNormEnd_xyz.c);
    free_dmatrix(pointNormEnd_yzb.m, 1, pointNormEnd_yzb.l, 1, pointNormEnd_yzb.c);
}

//Notice: I use not classical way to get normal vector
void sphere(unsigned int r, unsigned int g, unsigned int b)
{
    //information of sphere
    double radiumphere = 80;
    double v = 0.0 * M_PI;
    double u = 0.0 * M_PI;
    double dv = 2.0 * M_PI / 64;
    double du = M_PI / 32;
    double floatDepth_xyz, floatDepth_yzb, brightness_xyz, brightness_yzb;

    //initialization
    dmatrix_t x, y, z, buffer, center, vecYX, vecXZ, vecYZ, vecZbuffer, pointCenter_xyz, pointCenter_yzb, vecNorm_xyz, vecNorm_yzb, vecLight_xyz, vecLight_yzb, pointNormEnd_xyz, pointNormEnd_yzb;//x, y, z for torus, center of every circle on torus
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&vecYX, 4, 1);
    dmat_alloc(&vecXZ, 4, 1);
    dmat_alloc(&vecYZ, 4, 1);
    dmat_alloc(&vecZbuffer, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    dmat_alloc(&center, 4, 1);
    dmat_alloc(&pointCenter_xyz, 4, 1);
    dmat_alloc(&pointCenter_yzb, 4, 1);
    dmat_alloc(&vecNorm_xyz, 4, 1);
    dmat_alloc(&vecNorm_yzb, 4, 1);
    dmat_alloc(&vecLight_xyz, 4, 1);
    dmat_alloc(&vecLight_yzb, 4, 1);
    dmat_alloc(&pointNormEnd_xyz, 4, 1);
    dmat_alloc(&pointNormEnd_yzb, 4, 1);

    //the buffer is the x's opposite, so x, y, z, buffer makes a circle:
    // x - z
    // |
    // y - buffer
    buffer.m[1][1] = radiumphere * cos(v) * sin(u);
    buffer.m[2][1] = radiumphere * sin(v) * sin(u);
    buffer.m[3][1] = radiumphere * cos(u);
    buffer.m[4][1] = 1;

    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = -100.0;
    pointLight.m[2][1] = -100.0;
    pointLight.m[3][1] = 300.0;
    pointLight.m[4][1] = 1.0;
    for (v = 0 * M_PI; v - dv <= 2 * M_PI; v += dv)
    {
        for (u = 0.0; u - du <= M_PI; u += du)
        {
            //0.x, y, z calculation part;
            center.m[1][1] = 0;
            center.m[2][1] = 0;
            center.m[3][1] = 0;
            center.m[4][1] = 1;
            x.m[1][1] = radiumphere * cos(v) * sin(u);
            x.m[2][1] = radiumphere * sin(v) * sin(u);
            x.m[3][1] = radiumphere * cos(u);
            x.m[4][1] = 1;
            y.m[1][1] = radiumphere * cos(v) * sin(u + du);
            y.m[2][1] = radiumphere * sin(v) * sin(u + du);
            y.m[3][1] = radiumphere * cos(u + du);
            y.m[4][1] = 1;
            z.m[1][1] = radiumphere * cos(v + dv) * sin(u);
            z.m[2][1] = radiumphere * sin(v + dv) * sin(u);
            z.m[3][1] = radiumphere * cos(u);
            z.m[4][1] = 1;
            buffer.m[1][1] = radiumphere * cos(v + dv) * sin(u + du);
            buffer.m[2][1] = radiumphere * sin(v + dv) * sin(u + du);
            buffer.m[3][1] = radiumphere * cos(u + du);
            buffer.m[4][1] = 1;


            //1. edge vectors
            dmat_alloc(&vecYX, 4, 1);
            vecYX = *dmat_sub(&x, &y);//x-y
            vecXZ = *dmat_sub(&z, &x);
            vecYZ = *dmat_sub(&z, &y);
            vecZbuffer = *dmat_sub(&buffer, &z);


            //2.point center
            pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1]) / 3;
            pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1]) / 3;
            pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1]) / 3;
            pointCenter_xyz.m[4][1] = 1;

            pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1]) / 3;
            pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1]) / 3;
            pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1]) / 3;
            pointCenter_yzb.m[4][1] = 1;


            //3. Depth of each mesh
            floatDepth_xyz = dmat_norm(dmat_sub(&pointEye, &pointCenter_xyz));
            floatDepth_yzb = dmat_norm(dmat_sub(&pointEye, &pointCenter_yzb));


            //4.calculate Normal vectors
            //vecNorm_xyz = *crossProduct(vecYX, vecXZ); //the classical way
            vecNorm_xyz = *dmat_sub(&pointCenter_xyz, &center);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);
            vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.
            //vecNorm_yzb = *crossProduct(vecYZ, vecZbuffer); //the classical way
            vecNorm_yzb = *dmat_sub(&pointCenter_yzb, &center);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);
            vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.


            //5. light vector
            vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);


            //6. brightness: cos of two vector X light
            brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            brightness_yzb = brightness_xyz; //angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;


            //7.light-center & normal vector (optional)
            //7.1 xyz
            pointNormEnd_xyz = *dmat_add(&pointCenter_xyz,
                                         &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1],0,0);
            //normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1],0,0);

            //7.2 yzb
            pointNormEnd_yzb = *dmat_add(&pointCenter_yzb,
                                         &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1] ,0 ,0);
            //normal vector
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1],0,0);


            //8. draw the meshes.
            //Line3D(x, y, C);
            //Line3D(x, z, C);
            //Line3D(y, z, C);


            //9.fill with color
            x = *perspective_projection(dmat_mult(&C, &x));
            y = *perspective_projection(dmat_mult(&C, &y));
            z = *perspective_projection(dmat_mult(&C, &z));
            buffer = *perspective_projection(dmat_mult(&C, &buffer));

            XFillConvexPolygon(d, w, s, x, y, z, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz, floatDepth_xyz);
            XFillConvexPolygon(d, w, s, buffer, y, z, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb, floatDepth_yzb);
        }
    }


    //9.clean
    free_dmatrix(x.m, 1, x.l, 1, x.c);
    free_dmatrix(y.m, 1, y.l, 1, y.c);
    free_dmatrix(z.m, 1, z.l, 1, z.c);
    free_dmatrix(buffer.m, 1, buffer.l, 1, buffer.c);
    free_dmatrix(pointLight.m, 1, pointLight.l, 1, pointLight.c);
    free_dmatrix(vecYX.m, 1, vecYX.l, 1, vecYX.c);
    free_dmatrix(vecXZ.m, 1, vecXZ.l, 1, vecXZ.c);
    free_dmatrix(vecYZ.m, 1, vecYZ.l, 1, vecYZ.c);
    free_dmatrix(vecZbuffer.m, 1, vecZbuffer.l, 1, vecZbuffer.c);
    free_dmatrix(center.m, 1, center.l, 1, center.c);
    free_dmatrix(pointCenter_xyz.m, 1, pointCenter_xyz.l, 1, pointCenter_xyz.c);
    free_dmatrix(pointCenter_yzb.m, 1, pointCenter_yzb.l, 1, pointCenter_yzb.c);
    free_dmatrix(vecNorm_xyz.m, 1, vecNorm_xyz.l, 1, vecNorm_xyz.c);
    free_dmatrix(vecNorm_yzb.m, 1, vecNorm_yzb.l, 1, vecNorm_yzb.c);
    free_dmatrix(vecLight_xyz.m, 1, vecLight_xyz.l, 1, vecLight_xyz.c);
    free_dmatrix(vecLight_yzb.m, 1, vecLight_yzb.l, 1, vecLight_yzb.c);
    free_dmatrix(pointNormEnd_xyz.m, 1, pointNormEnd_xyz.l, 1, pointNormEnd_xyz.c);
    free_dmatrix(pointNormEnd_yzb.m, 1, pointNormEnd_yzb.l, 1, pointNormEnd_yzb.c);

}

//void cone(unsigned int red, unsigned int g, unsigned int b)
//{
//    //matrix designed for cone;
//    // dmatrix_t D = C;
//
//    // D = *rotation('z', M_PI/12, D);
//    // //D = *rotation('x', 12*M_PI/12, D);
//    // //D = *rotation('y', -1*M_PI/12, D);
//    // D = *translation(0, 0, -480, D);
//
//
//    //information of sphere
//    double floatDepth_xyz, floatDepth_yzb, brightness_xyz, brightness_yzb;
//    double height = 200;
//    double r = 10;
//    double theta = 0.0 * M_PI;     //range of sphere
//    double u = 0.0;
//    double dtheta = M_PI / 16;
//    double du = height / 10;
//
//    //initialization
//    dmatrix_t P[3];
//    for (int i = 0; i < 3; i++) dmat_alloc(&P[i], 4, 1);
//    dmatrix_t x, y, z, buffer, center, vecYX, vecXZ, vecYZ, vecZbuffer, pointCenter_xyz, pointCenter_yzb, vecNorm_xyz, vecNorm_yzb, vecLight_xyz, vecLight_yzb, pointNormEnd_xyz, pointNormEnd_yzb;//x, y, z for torus, center of every circle on torus
//    dmat_alloc(&x, 4, 1);
//    dmat_alloc(&y, 4, 1);
//    dmat_alloc(&z, 4, 1);
//    dmat_alloc(&vecYX, 4, 1);
//    dmat_alloc(&vecXZ, 4, 1);
//    dmat_alloc(&vecYZ, 4, 1);
//    dmat_alloc(&vecZbuffer, 4, 1);
//    dmat_alloc(&buffer, 4, 1);
//    dmat_alloc(&center, 4, 1);
//    dmat_alloc(&pointCenter_xyz, 4, 1);
//    dmat_alloc(&pointCenter_yzb, 4, 1);
//    dmat_alloc(&vecNorm_xyz, 4, 1);
//    dmat_alloc(&vecNorm_yzb, 4, 1);
//    dmat_alloc(&vecLight_xyz, 4, 1);
//    dmat_alloc(&vecLight_yzb, 4, 1);
//    dmat_alloc(&pointNormEnd_xyz, 4, 1);
//    dmat_alloc(&pointNormEnd_yzb, 4, 1);
//
//    //the buffer is the x's opposite, so x, y, z, buffer makes a circle:
//    // x - z
//    // |
//    // y - buffer
//    buffer.m[1][1] = r * cos(theta) * (height - u) / height;
//    buffer.m[2][1] =  r * sin(theta) * (height - u) / height;
//    buffer.m[3][1] = u;
//    buffer.m[4][1] = 1;
//    buffer = *perspective_projection(dmat_mult(&C, &buffer));
//
//    //light
//    dmatrix_t pointLight;
//    dmat_alloc(&pointLight, 4, 1);
//    pointLight.m[1][1] = -100.0;
//    pointLight.m[2][1] = -100.0;
//    pointLight.m[3][1] = 300.0;
//    pointLight.m[4][1] = 1.0;
//
//    /* Notice: the sphere's u, that is every slices' equation, only need run a half */
//    double delta = 2; double start = 0;
//    for (u = 0; u + du <= height; u += du)
//    {
//        for (theta = (start) * M_PI; theta <= (delta + start) * M_PI; theta += dtheta)
//        {
//            //0.x, y, z calculation part;
//            x.m[1][1] = r * cos(theta) * (u) / height;
//            x.m[2][1] = r * sin(theta) * (u) / height;
//            x.m[3][1] = u;
//            x.m[4][1] = 1;
//            y.m[1][1] = r * cos(theta + dtheta) * (u) / height;
//            y.m[2][1] = r * sin(theta + dtheta) * (u) / height;
//            y.m[3][1] = u;
//            y.m[4][1] = 1;
//            z.m[1][1] = r * cos(theta) * (u + du) / height;
//            z.m[2][1] = r * sin(theta) * (u + du) / height;
//            z.m[3][1] = u + du;
//            z.m[4][1] = 1;
//            buffer.m[1][1] = r * cos(theta + dtheta) * (u + du) / height;
//            buffer.m[2][1] = r * sin(theta + dtheta) *  (u + du) / height;
//            buffer.m[3][1] = u + du;
//            buffer.m[4][1] = 1;
//
//            //1. edge vectors
//            dmat_alloc(&vecYX, 4, 1);
//            vecYX = *dmat_sub(&x, &y);//x-y
//            vecXZ = *dmat_sub(&z, &x);
//            vecYZ = *dmat_sub(&z, &y);
//            vecZbuffer = *dmat_sub(&buffer, &z);
//
//            //2.point center
//            pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1]) / 3;
//            pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1]) / 3;
//            pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1]) / 3;
//            pointCenter_xyz.m[4][1] = 1;
//
//            pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1]) / 3;
//            pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1]) / 3;
//            pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1]) / 3;
//            pointCenter_yzb.m[4][1] = 1;
//
//
//            //3. Depth of each mesh
//            floatDepth_xyz = dmat_norm(dmat_sub(&pointEye, &pointCenter_xyz));
//            floatDepth_yzb = dmat_norm(dmat_sub(&pointEye, &pointCenter_yzb));
//
//
//            //4.calculate Normal vectors
//            //vecNorm_xyz = *crossProduct(vecYX, vecXZ); //the classical way
//            vecNorm_xyz = *dmat_sub(&pointCenter_xyz, &center);
//            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);
//            vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.
//
//            //vecNorm_yzb = *crossProduct(vecYZ, vecZbuffer); //the classical way
//            vecNorm_yzb = *dmat_sub(&pointCenter_yzb, &center);
//            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);
//            vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.
//
//
//            //5. light vector
//            vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
//            vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);
//
//
//            //6. brightness: cos of two vector X light
//            brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
//            brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
//
//
//            //7.light-center & normal vector (optional)
//            //7.1 xyz
//            pointNormEnd_xyz = *dmat_add(&pointCenter_xyz,
//                                         &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
//            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
//            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
//            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
//            //light to point
//            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1], 0, 0);
//            //normal vector
//            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1],0,0);
//
//            //7.2 yzb
//            pointNormEnd_yzb = *dmat_add(&pointCenter_yzb,
//                                         &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
//            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
//            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
//            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
//            //light to point
//            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1], 0, 0);
//            //normal vector
//            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1],0,0);
//
//
//            //8. draw the meshes.
//            //Line3D(x, y, C);
//            //Line3D(x, z, C);
//            //Line3D(y, z, C);
//
//
//            //9.fill with color
//            x = *perspective_projection(dmat_mult(&C, &x));
//            y = *perspective_projection(dmat_mult(&C, &y));
//            z = *perspective_projection(dmat_mult(&C, &z));
//            buffer = *perspective_projection(dmat_mult(&C, &buffer));
//            P[0] = x;
//            P[1] = y;
//            P[2] = z;
//            XFillConvexPolygon(d, w, s, P, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz, floatDepth_xyz);
//            P[0] = buffer;
//            XFillConvexPolygon(d, w, s, P, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb, floatDepth_yzb);
//        }
//    }
//    //9.clean
//    free_dmatrix(x.m, 1, x.l, 1, x.c);
//    free_dmatrix(y.m, 1, y.l, 1, y.c);
//    free_dmatrix(z.m, 1, z.l, 1, z.c);
//    free_dmatrix(buffer.m, 1, buffer.l, 1, buffer.c);
//    free_dmatrix(pointLight.m, 1, pointLight.l, 1, pointLight.c);
//    free_dmatrix(vecYX.m, 1, vecYX.l, 1, vecYX.c);
//    free_dmatrix(vecXZ.m, 1, vecXZ.l, 1, vecXZ.c);
//    free_dmatrix(vecYZ.m, 1, vecYZ.l, 1, vecYZ.c);
//    free_dmatrix(vecZbuffer.m, 1, vecZbuffer.l, 1, vecZbuffer.c);
//    free_dmatrix(center.m, 1, center.l, 1, center.c);
//    free_dmatrix(pointCenter_xyz.m, 1, pointCenter_xyz.l, 1, pointCenter_xyz.c);
//    free_dmatrix(pointCenter_yzb.m, 1, pointCenter_yzb.l, 1, pointCenter_yzb.c);
//    free_dmatrix(vecNorm_xyz.m, 1, vecNorm_xyz.l, 1, vecNorm_xyz.c);
//    free_dmatrix(vecNorm_yzb.m, 1, vecNorm_yzb.l, 1, vecNorm_yzb.c);
//    free_dmatrix(vecLight_xyz.m, 1, vecLight_xyz.l, 1, vecLight_xyz.c);
//    free_dmatrix(vecLight_yzb.m, 1, vecLight_yzb.l, 1, vecLight_yzb.c);
//    free_dmatrix(pointNormEnd_xyz.m, 1, pointNormEnd_xyz.l, 1, pointNormEnd_xyz.c);
//    free_dmatrix(pointNormEnd_yzb.m, 1, pointNormEnd_yzb.l, 1, pointNormEnd_yzb.c);
//}

//
//void meteor(unsigned int r, unsigned int g, unsigned int b)
//{
//    //information of sphere
//    double radiumphere = 80;
//
//    //range of sphere
//    double v = 0.0 * M_PI;
//    double u = 0.0 * M_PI;
//    double dv = 2.0 * M_PI / 80;
//    double du = M_PI / 32;
//
//    //light
//    dmatrix_t pointLight;
//    dmat_alloc(&pointLight, 4, 1);
//    pointLight.m[1][1] = 180.0;
//    pointLight.m[2][1] = 180.0;
//    pointLight.m[3][1] = 180.0;
//    pointLight.m[4][1] = 1.0;
//
//    //points to draw sphere. the buffer here has same idea as in Torus'
//    dmatrix_t P[3];
//    for (int i = 0; i < 3; i++)
//    {
//        dmat_alloc(&P[i], 4, 1);
//    }
//    dmatrix_t x, y, z, buffer;
//    dmat_alloc(&x, 4, 1);
//    dmat_alloc(&y, 4, 1);
//    dmat_alloc(&z, 4, 1);
//    dmat_alloc(&buffer, 4, 1);
//    buffer.m[1][1] = radiumphere * cos(v) * sin(u);
//    buffer.m[2][1] = radiumphere * sin(v) * sin(u);
//    buffer.m[3][1] = radiumphere * cos(u);
//    buffer.m[4][1] = 1;
//    buffer = *perspective_projection(dmat_mult(&C, &buffer));
//
//
//    /* Notice: the sphere's u, that is every slices' euqation, only need run a half */
//    double delta = 2;
//    double start = 0;
//    for (v = (start) * M_PI; v + dv <= (delta + start) * M_PI; v += dv)
//    {
//        for (u = 0.0; u + du <= M_PI; u += du)
//        {
//            x.m[1][1] = radiumphere * cos(v) * sin(u);
//            x.m[2][1] = radiumphere * sin(v) * sin(u);
//            x.m[3][1] = radiumphere * cos(u);
//            x.m[4][1] = 1;
//            y.m[1][1] = radiumphere * cos(v) * sin(u + du);
//            y.m[2][1] = radiumphere * sin(v) * sin(u + du);
//            y.m[3][1] = radiumphere * cos(u + du);
//            y.m[4][1] = 1;
//            z.m[1][1] = radiumphere * cos(v + dv) * sin(u);
//            z.m[2][1] = radiumphere * sin(v + dv) * sin(u);
//            z.m[3][1] = radiumphere * sin(u);
//            z.m[4][1] = 1;
//            buffer.m[1][1] = radiumphere * cos(v + dv) * sin(u + du);
//            buffer.m[2][1] = radiumphere * sin(v + dv) * sin(u + du);
//            buffer.m[3][1] = radiumphere * cos(u + du);
//            buffer.m[4][1] = 1;
//
//            //2.light and calculation of brightness;
//            //calcuate vectors based on surface x-y-z-buffer
//            dmatrix_t vecYX;
//            dmat_alloc(&vecYX, 4, 1);
//            vecYX = *dmat_sub(&x, &y);//x-y
//            dmatrix_t vecXZ;
//            dmat_alloc(&vecXZ, 4, 1);
//            vecXZ = *dmat_sub(&z, &x);
//            dmatrix_t vecYZ;
//            dmat_alloc(&vecYZ, 4, 1);
//            vecYZ = *dmat_sub(&z, &y);
//            dmatrix_t vecZbuffer;
//            dmat_alloc(&vecZbuffer, 4, 1);
//            vecZbuffer = *dmat_sub(&buffer, &z);
//
//            //calcuate two vectors: xyz, yzbuffer
//            dmatrix_t vecNorm_xyz;
//            dmat_alloc(&vecNorm_xyz, 4, 1);
//            vecNorm_xyz = *crossProduct(vecXZ, vecYX);
//            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);
//            vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz,
//                                            100); //normailize it and times 10 could make it looks clear when visualize normalvector.
//
//            dmatrix_t vecNorm_yzb;
//            dmat_alloc(&vecNorm_yzb, 4, 1);
//            vecNorm_yzb = *crossProduct(vecZbuffer, vecYZ);
//            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);
//            vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb,
//                                            100); //normailize it and times 10 could make it looks clear when visualize normalvector.
//
//            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
//            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);
//
//            //calcuate center of xyz, yzb
//            dmatrix_t pointCenter_xyz;
//            dmat_alloc(&pointCenter_xyz, 4, 1);
//            pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1]) / 3;
//            pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1]) / 3;
//            pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1]) / 3;
//            pointCenter_xyz.m[4][1] = 1;
//
//            dmatrix_t pointCenter_yzb;
//            dmat_alloc(&pointCenter_yzb, 4, 1);
//            pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1]) / 3;
//            pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1]) / 3;
//            pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1]) / 3;
//            pointCenter_yzb.m[4][1] = 1;
//
//            //calculate two light vector
//            dmatrix_t vecLight_xyz;
//            dmat_alloc(&vecLight_xyz, 4, 1);
//            vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
//            dmatrix_t vecLight_yzb;
//            dmat_alloc(&vecLight_yzb, 4, 1);
//            vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);
//
//            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
//            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
//            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
//            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);
//
//            //calculate cos of two vector X light
//            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
//            double brightness_xyz = angle(vecNorm_xyz, vecLight_xyz);
//            if (brightness_xyz < 0)brightness_xyz = -brightness_xyz;
//            //double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
//            double brightness_yzb = angle(vecNorm_yzb, vecLight_yzb);
//            if (brightness_yzb < 0)brightness_yzb = -brightness_yzb;
//
//
//            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;
//
//
//            //3.draw centre to light part and normal vector(optional)
//            //xyz
//            dmatrix_t pointNormEnd_xyz;
//            dmat_alloc(&pointNormEnd_xyz, 4, 1);
//            pointNormEnd_xyz = *dmat_add(&pointCenter_xyz,
//                                         &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
//            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
//            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
//            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
//            //light to point, and normal vector
//            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
//            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);
//
//            //yzb
//            dmatrix_t pointNormEnd_yzb;
//            dmat_alloc(&pointNormEnd_yzb, 4, 1);
//            pointNormEnd_yzb = *dmat_add(&pointCenter_yzb,
//                                         &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
//            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
//            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
//            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
//            //light to point, and normal vector
//            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
//            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1]);
//
//
//            //4. draw the meshes.
//            x = *perspective_projection(dmat_mult(&C, &x));
//            y = *perspective_projection(dmat_mult(&C, &y));
//            z = *perspective_projection(dmat_mult(&C, &z));
//            buffer = *perspective_projection(dmat_mult(&C, &buffer));
//
//            //Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
//            //Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
//            //Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
//            //5.fill with color
//            P[0] = x;
//            P[1] = y;
//            P[2] = z;
//            XFillConvexPolygon(d, w, s, P, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz);
//            P[0] = buffer;
//            XFillConvexPolygon(d, w, s, P, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb);
//        }
//    }
//    //7.clean
//    free_dmatrix(x.m, 1, x.l, 1, x.c);
//    free_dmatrix(y.m, 1, y.l, 1, y.c);
//    free_dmatrix(z.m, 1, z.l, 1, z.c);
//    free_dmatrix(buffer.m, 1, buffer.l, 1, buffer.c);
//    free_dmatrix(pointLight.m, 1, pointLight.l, 1, pointLight.c);
//}
//
//
//void plane_cut_and_fill(dmatrix_t x, dmatrix_t y, dmatrix_t z, dmatrix_t D, unsigned int r, unsigned int g, unsigned int b )
//{
//    if (dmat_compare(&x, &z)==true)
//        return;
//
//    dmatrix_t dxy, dxz, dyz, u, v, v1, v2, pointCenter_xyz, pointCenter_yzb;//x, y, z for plane
//    u = *dmat_duplicate(&x);
//
//    double delta  = 0.1;
//    dxy = *dmat_scalar_mult(dmat_sub(&y,&x),delta);
//    dxz = *dmat_scalar_mult(dmat_sub(&z,&x),delta);
//    dyz = *dmat_scalar_mult(dmat_sub(&z,&y),delta);
//
////        X
////        u
////       / \
////      v - v1
////       \/
////   Y   v2   Z
//
//    for(u; dmat_compare(&y, dmat_add(&u, &dxy))==false; dmat_add(&u, &dxy))//这样应该会缺左下一个角
//    {
//        v = *dmat_add(&u, &dxy);
//        v1 = *dmat_add(&u, &dxy);
//        v2 = *dmat_add(&v1, &dxy);
//
//        pointCenter_xyz = *dmat_add(&u, &v);
//        pointCenter_xyz = *dmat_add(&pointCenter_xyz, &v1);
//        pointCenter_xyz = *dmat_scalar_mult(&pointCenter_xyz, 1/3);
//        double floatDepth_xyz = dmat_norm(dmat_sub(&pointEye, &pointCenter_xyz));
//
//
//        pointCenter_yzb = *dmat_add(&v, &v1);
//        pointCenter_yzb = *dmat_add(&pointCenter_yzb, &v2);
//        pointCenter_yzb = *dmat_scalar_mult(&pointCenter_yzb, 1/3);
//        double floatDepth_yzb = dmat_norm(dmat_sub(&pointEye, &pointCenter_yzb));
//
//        Line3D(u, v, C);
//        Line3D(v2, v, C);
//        Line3D(v1, v, C);  //So the right most edge is missed;
//
//        dmatrix_t P[3];
//        for (int i = 0; i < 3; i++)dmat_alloc(&P[i], 4, 1);
//        u = *perspective_projection(dmat_mult(&D, &u));
//        v = *perspective_projection(dmat_mult(&D, &v));
//        v1 = *perspective_projection(dmat_mult(&D, &v1));
//        v2 = *perspective_projection(dmat_mult(&D, &v2));
//        P[0] = u;
//        P[1] = v;
//        P[2] = v1;
//        XFillConvexPolygon(d, w, s, P, 3, r, g, b, floatDepth_xyz);
//        P[0] = v2;
//        XFillConvexPolygon(d, w, s, P, 3, r, g, b, floatDepth_yzb);
//    }
//
//    plane_cut_and_fill(*dmat_add(&x, &dxz), *dmat_add(&y, &dyz), z, D, r, g, b);
//
//
//}
//
//void plane_write_up()
//{
//    dmatrix_t x, y, z, pointCenter_xyz;//x, y, z for plane
//    dmat_alloc(&x, 4, 1);
//    dmat_alloc(&y, 4, 1);
//    dmat_alloc(&z, 4, 1);
//    dmat_alloc(&pointCenter_xyz, 4, 1);
//
//    x.m[1][1] = 0;
//    x.m[2][1] = 0;
//    x.m[3][1] = 0;
//    x.m[4][1] = 1;
//
//    y.m[1][1] = 200;
//    y.m[2][1] = 200;
//    y.m[3][1] = 200;
//    y.m[4][1] = 1;
//
//    z.m[1][1] = 200;
//    z.m[2][1] = 0;
//    z.m[3][1] = 100;
//    z.m[4][1] = 1;
//
//    pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1]) / 3;
//    pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1]) / 3;
//    pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1]) / 3;
//    pointCenter_xyz.m[4][1] = 1;
//    double floatDepth_xyz = dmat_norm(dmat_sub(&pointEye, &pointCenter_xyz));
//    cout<<"plane write up's distance: "<<floatDepth_xyz<<endl;
//
//    Line3D(x, y, C);
//    Line3D(x, z, C);
//    Line3D(z, y, C);
//    dmatrix_t P[3];
//    for (int i = 0; i < 3; i++)dmat_alloc(&P[i], 4, 1);
//    x = *perspective_projection(dmat_mult(&C, &x));
//    y = *perspective_projection(dmat_mult(&C, &y));
//    z = *perspective_projection(dmat_mult(&C, &z));
//    P[0] = x;
//    P[1] = y;
//    P[2] = z;
//
//
//    XFillConvexPolygon(d, w, s, P, 3, 250, 250, 250, floatDepth_xyz);
//}
//
//void plane_dark_down()
//{
//    dmatrix_t x, y, z, dxy, dxz, u, v, pointCenter_xyz;//x, y, z for plane
//    dmat_alloc(&x, 4, 1);
//    dmat_alloc(&y, 4, 1);
//    dmat_alloc(&z, 4, 1);
//    dmat_alloc(&u, 4, 1);
//    dmat_alloc(&v, 4, 1);
//    dmat_alloc(&dxy, 4, 1);
//    dmat_alloc(&dxz, 4, 1);
//    dmat_alloc(&pointCenter_xyz, 4, 1);
//
//    x.m[1][1] = 0;
//    x.m[2][1] = 0;
//    x.m[3][1] = -200;
//    x.m[4][1] = 1;
//
//    y.m[1][1] = 100;
//    y.m[2][1] = 100;
//    y.m[3][1] = 50;
//    y.m[4][1] = 1;
//
//    z.m[1][1] = 200;
//    z.m[2][1] = 0;
//    z.m[3][1] = 50;
//    z.m[4][1] = 1;
//
//    plane_cut_and_fill(x, y, z, C, 100, 100, 100);
//}

void x_axi (dmatrix_t camera)
{
    dmatrix_t o, x;
    dmat_alloc(&o, 4, 1);
    dmat_alloc(&x, 4, 1);

    o.m[1][1] = 0;
    o.m[2][1] = 0;
    o.m[3][1] = 0;
    o.m[4][1] = 1;
    x.m[1][1] = 200;
    x.m[2][1] = 0;
    x.m[3][1] = 0;
    x.m[4][1] = 1;

    Line3D(o, x, camera);
}

void y_axi(dmatrix_t camera)
{
    dmatrix_t o, y;
    dmat_alloc(&o, 4, 1);
    dmat_alloc(&y, 4, 1);

    o.m[1][1] = 0;
    o.m[2][1] = 0;
    o.m[3][1] = 0;
    o.m[4][1] = 1;
    y.m[1][1] = 0;
    y.m[2][1] = 200;
    y.m[3][1] = 0;
    y.m[4][1] = 1;

    Line3D(o, y, camera);
}

void z_axi(dmatrix_t camera)
{
    dmatrix_t o, z;//x, y, z for plane
    dmat_alloc(&o, 4, 1);
    dmat_alloc(&z, 4, 1);

    o.m[1][1] = 0;
    o.m[2][1] = 0;
    o.m[3][1] = 0;
    o.m[4][1] = 1;
    z.m[1][1] = 0;
    z.m[2][1] = 0;
    z.m[3][1] = -200;
    z.m[4][1] = 1;

    Line3D(o, z, camera);
}

void original(dmatrix_t camera)
{
    dmatrix_t o;//x, y, z for plane
    dmat_alloc(&o, 4, 1);
    o.m[1][1] = 0;
    o.m[2][1] = 0;
    o.m[3][1] = 0;
    o.m[4][1] = 1;
    o = *perspective_projection(dmat_mult(&camera, &o));
    DrawPixel(o.m[1][1], o.m[2][1], 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1] + 1, o.m[2][1], 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1] - 1, o.m[2][1], 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1], o.m[2][1] + 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1], o.m[2][1] - 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1] + 1, o.m[2][1] + 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1] - 1, o.m[2][1] - 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1] + 1, o.m[2][1] - 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(o.m[1][1] - 1, o.m[2][1] + 1, 255, 0, 0, 0);//original is always the closest;
}

void coordinate(bool x, bool y, bool z, bool origin)
{
    if (x) x_axi(C_original);
    if (y) y_axi(C_original);
    if (z) z_axi(C_original);
    if (origin) original(C_original);
}

void eyes()
{
    dmatrix_t o, eye;//x, y, z for plane
    eye = pointEye;

    dmat_alloc(&o, 4, 1);
    o.m[1][1] = 0;
    o.m[2][1] = 0;
    o.m[3][1] = 0;
    o.m[4][1] = 1;

    Line3D(o, eye, C);

    //eye = *perspective_projection(dmat_mult(&C, &eye));
    DrawPixel(eye.m[1][1], eye.m[2][1], 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1] + 1, eye.m[2][1], 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1] - 1, eye.m[2][1], 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1], eye.m[2][1] + 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1], eye.m[2][1] - 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1] + 1, eye.m[2][1] + 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1] - 1, eye.m[2][1] - 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1] + 1, eye.m[2][1] - 1, 255, 0, 0, 0);//original is always the closest;
    DrawPixel(eye.m[1][1] - 1, eye.m[2][1] + 1, 255, 0, 0, 0);//original is always the closest;
}

void Draw()
{
    //cone(200, 200, 200);
    sphere(237, 189, 101);
    torus(220, 220, 220);
    //plane_write_up();
    //plane_dark_down();

    //eyes();
    //coordinate(true, true, true, true);
}

void OnDisplay()
{
    for (long i = 0; i < windowH * windowW*3; i++){ frame.push_back(255); }
    for (long i = 0; i < windowH * windowW; i++) { depth.push_back(99999999); }
    //main print function
    Draw();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(windowW, windowH, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte *)&frame[0]);
    glFlush();
    glutSwapBuffers();

    vector<char>().swap(frame);
    vector<double>().swap(depth);
}


