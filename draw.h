
#include <X11/Xlib.h>
#include <OpenGL/gl.h>
#include <OpenGl/glu.h>
#include <GLUT/glut.h>
#include <string.h>
#include "camera.h"

int REMOVAL = 1;
const int X = 1;
const int Y = 2;

double angle_total = 0.0;  //totally rotated angle

Display *d;
Window w;
XEvent e;
int s;

unsigned char frame[windowW * windowH * 3];
unsigned char window[windowW * windowH];

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

void hiddenremover()
{
    memset(window, 255, windowW * windowH);
    REMOVAL = 0;
}

void OnKeyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 'q':
            printf("total rotate angle: %f\n",angle_total);
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

    if (0 <= dy && dy <= dx)
    { // 0 < m < 1, case 2
        //std::cout << "case 2" << std::endl;
        for (int x = x1; x <= x2; x += 1)
        {
            DrawPixel(x, y, r, g, b);
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
        for (int x = x1; x >= x2; x -= 1)
        {
            DrawPixel(x, y, r, g, b);
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
        for (int x = x1; x <= x2; x += 1)
        {
            DrawPixel(x, y, r, g, b);
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
        Bresenham(x2, y2, x1, y1, r, g, b);
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

        int x = x1;
        for (int y = y1; y <= y2; y += 1)
        {
            DrawPixel(x, y, r, g, b);
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

        dx = -dx;
        int x = x2;
        for (int y = y2; y >= y1; y -= 1)
        {
            DrawPixel(x, y, r, g, b);
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
        Bresenham(x2, y2, x1, y1, r, g, b);
    } else if (0 <= dx && dy <= -dx)
    { // m < -1, case 4
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
    Bresenham((int) x1, (int) y1, (int) x2, (int) y2, 0, 0, 0);
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
    window[1] = 0;
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

void XFillConvexPolygon(Display *d, Window w, int s, dmatrix_t P[], int n, unsigned int r, unsigned int g, unsigned int bb)
{

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
            DrawPixel(i, y, r, g, bb);
        }
    }
    free(horizontal);
    free(active);
    free(intersections);
}

void torus(unsigned int r, unsigned int g, unsigned int b)
{
    //information of torus
    double radium_torus, radium_torus_pipe;
    radium_torus = 180;
    radium_torus_pipe = 2;
    //the range of torus
    double v = 0.0 * M_PI;
    double u = 0.0 * M_PI;
    double du = 2.0 * M_PI / 256;
    double dv = 2.0 * M_PI / 256;

    //points to draw the torus
    dmatrix_t P[3];
    for (int i = 0; i < 3; i++)
    {
        dmat_alloc(&P[i], 4, 1);
    }
    dmatrix_t x, y, z, buffer;//x, y, z for torus
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    //the buffer is the x's opposite, so x, y, z, buffer makes a circle:
    // x - z
    // |
    // y - buffer
    buffer.m[1][1] = (radium_torus + radium_torus_pipe * cos(u)) * cos(v);
    buffer.m[2][1] = (radium_torus + radium_torus_pipe * cos(u)) * sin(v);
    buffer.m[3][1] = radium_torus_pipe * sin(u);;
    buffer.m[4][1] = 1;
    buffer = *perspective_projection(dmat_mult(&C, &buffer));

    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = 180.0;
    pointLight.m[2][1] = 180.0;
    pointLight.m[3][1] = 180.0;
    pointLight.m[4][1] = 1.0;

    for (v = 0 * M_PI; v <= 2.0 * M_PI; v += dv)
    {
        for (u = 0 * M_PI; u <= 1.0 * M_PI; u += du)
        {
    //1.x, y, z calculation part;
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

    //2.light and calculation of brightness;
            //calcuate vectors based on surface x-y-z-buffer
            dmatrix_t vecYX;dmat_alloc(&vecYX, 4, 1);vecYX = *dmat_sub(&x, &y);//x-y
            dmatrix_t vecXZ;dmat_alloc(&vecXZ, 4, 1);vecXZ = *dmat_sub(&z, &x);
            dmatrix_t vecYZ;dmat_alloc(&vecYZ, 4, 1);vecYZ = *dmat_sub(&z, &y);
            dmatrix_t vecZbuffer;dmat_alloc(&vecZbuffer, 4, 1);vecZbuffer = *dmat_sub(&buffer, &z);

            //calcuate two vectors: xyz, yzbuffer
            dmatrix_t vecNorm_xyz; dmat_alloc(&vecNorm_xyz, 4, 1);vecNorm_xyz = *crossProduct(vecYX, vecXZ);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            dmatrix_t vecNorm_yzb; dmat_alloc(&vecNorm_yzb, 4, 1);vecNorm_yzb = *crossProduct(vecYZ, vecZbuffer);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);

            //calcuate center of xyz, yzb
            dmatrix_t pointCenter_xyz;dmat_alloc(&pointCenter_xyz,4,1);pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1])/3;pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1])/3;pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1])/3;pointCenter_xyz.m[4][1] = 1;

            dmatrix_t pointCenter_yzb;dmat_alloc(&pointCenter_yzb,4,1);pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1])/3;pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1])/3;pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1])/3;pointCenter_yzb.m[4][1] = 1;

            //calculate two light vector
            dmatrix_t vecLight_xyz;dmat_alloc(&vecLight_xyz, 4, 1);vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            dmatrix_t vecLight_yzb;dmat_alloc(&vecLight_yzb, 4, 1);vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);

            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);

            //calculate cos of two vector X light
            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = -brightness_xyz;
            //double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
            double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = -brightness_yzb;

            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;


    //3.draw centre to light part and normal vector(optional)
            //xyz
            dmatrix_t pointNormEnd_xyz; dmat_alloc(&pointNormEnd_xyz, 4, 1);pointNormEnd_xyz = *dmat_add(&pointCenter_xyz, &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);

            //yzb
            dmatrix_t pointNormEnd_yzb; dmat_alloc(&pointNormEnd_yzb, 4, 1);pointNormEnd_yzb = *dmat_add(&pointCenter_yzb, &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1]);


    //4. draw the meshes.
            x = *perspective_projection(dmat_mult(&C, &x));
            y = *perspective_projection(dmat_mult(&C, &y));
            z = *perspective_projection(dmat_mult(&C, &z));
            buffer = *perspective_projection(dmat_mult(&C, &buffer));

            //Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
            //Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
            //Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
    //5.fill with color
            P[0] = x;P[1] = y;P[2] = z;
                        XFillConvexPolygon(d, w, s, P, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz);
            P[0] = buffer;
                        XFillConvexPolygon(d, w, s, P, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb);
        }
    }

    //6.clean
    free_dmatrix(x.m,1,x.l,1,x.c) ;
    free_dmatrix(y.m,1,y.l,1,y.c) ;
    free_dmatrix(z.m,1,z.l,1,z.c) ;
    free_dmatrix(buffer.m,1,buffer.l,1,buffer.c) ;
    free_dmatrix(pointLight.m,1,pointLight.l,1,pointLight.c) ;
}
void sphere(unsigned int r, unsigned int g, unsigned int b)
{
    //information of sphere
    double radiumphere = 80;

    //range of sphere
    double v = 0.0 * M_PI;
    double u = 0.0 * M_PI;
    double dv = 2.0 * M_PI / 80;
    double du = M_PI / 32;

    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = 180.0;
    pointLight.m[2][1] = 180.0;
    pointLight.m[3][1] = 180.0;
    pointLight.m[4][1] = 1.0;

    //points to draw sphere. the buffer here has same idea as in Torus'
    dmatrix_t P[3];
    for (int i = 0; i < 3; i++)
    {
        dmat_alloc(&P[i], 4, 1);
    }
    dmatrix_t x, y, z, buffer;
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    buffer.m[1][1] = radiumphere * cos(v) * sin(u);
    buffer.m[2][1] = radiumphere * sin(v) * sin(u);
    buffer.m[3][1] = radiumphere * cos(u);
    buffer.m[4][1] = 1;
    buffer = *perspective_projection(dmat_mult(&C, &buffer));


    /* Notice: the sphere's u, that is every slices' euqation, only need run a half */
    double delta = 2;
    double start = 0;
    for (v = (start) * M_PI; v+dv <= (delta+start) * M_PI; v += dv)
    {
        for (u = 0.0; u+du <= M_PI; u += du)
        {
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

    //2.light and calculation of brightness;
            //calcuate vectors based on surface x-y-z-buffer
            dmatrix_t vecYX;dmat_alloc(&vecYX, 4, 1);vecYX = *dmat_sub(&x, &y);//x-y
            dmatrix_t vecXZ;dmat_alloc(&vecXZ, 4, 1);vecXZ = *dmat_sub(&z, &x);
            dmatrix_t vecYZ;dmat_alloc(&vecYZ, 4, 1);vecYZ = *dmat_sub(&z, &y);
            dmatrix_t vecZbuffer;dmat_alloc(&vecZbuffer, 4, 1);vecZbuffer = *dmat_sub(&buffer, &z);

            //calcuate two vectors: xyz, yzbuffer
            dmatrix_t vecNorm_xyz; dmat_alloc(&vecNorm_xyz, 4, 1);vecNorm_xyz = *crossProduct(vecXZ, vecYX);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            dmatrix_t vecNorm_yzb; dmat_alloc(&vecNorm_yzb, 4, 1);vecNorm_yzb = *crossProduct(vecZbuffer, vecYZ);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);

            //calcuate center of xyz, yzb
            dmatrix_t pointCenter_xyz;dmat_alloc(&pointCenter_xyz,4,1);pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1])/3;pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1])/3;pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1])/3;pointCenter_xyz.m[4][1] = 1;

            dmatrix_t pointCenter_yzb;dmat_alloc(&pointCenter_yzb,4,1);pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1])/3;pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1])/3;pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1])/3;pointCenter_yzb.m[4][1] = 1;

            //calculate two light vector
            dmatrix_t vecLight_xyz;dmat_alloc(&vecLight_xyz, 4, 1);vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            dmatrix_t vecLight_yzb;dmat_alloc(&vecLight_yzb, 4, 1);vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);

            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);

            //calculate cos of two vector X light
            double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = -brightness_xyz;
            double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
            //double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = -brightness_yzb;


            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;


    //3.draw centre to light part and normal vector(optional)
            //xyz
            dmatrix_t pointNormEnd_xyz; dmat_alloc(&pointNormEnd_xyz, 4, 1);pointNormEnd_xyz = *dmat_add(&pointCenter_xyz, &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);

            //yzb
            dmatrix_t pointNormEnd_yzb; dmat_alloc(&pointNormEnd_yzb, 4, 1);pointNormEnd_yzb = *dmat_add(&pointCenter_yzb, &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1]);


    //4. draw the meshes.
            x = *perspective_projection(dmat_mult(&C, &x));
            y = *perspective_projection(dmat_mult(&C, &y));
            z = *perspective_projection(dmat_mult(&C, &z));
            buffer = *perspective_projection(dmat_mult(&C, &buffer));

            //Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
            //Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
            //Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
    //5.fill with color
            P[0] = x;P[1] = y;P[2] = z;
            XFillConvexPolygon(d, w, s, P, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz);
            P[0] = buffer;
            XFillConvexPolygon(d, w, s, P, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb);
        }
    }
    //7.clean
    free_dmatrix(x.m,1,x.l,1,x.c) ;
    free_dmatrix(y.m,1,y.l,1,y.c) ;
    free_dmatrix(z.m,1,z.l,1,z.c) ;
    free_dmatrix(buffer.m,1,buffer.l,1,buffer.c) ;
    free_dmatrix(pointLight.m,1,pointLight.l,1,pointLight.c) ;
}
void meteor(unsigned int r, unsigned int g, unsigned int b)
{
    //information of sphere
    double radiumphere = 80;

    //range of sphere
    double v = 0.0 * M_PI;
    double u = 0.0 * M_PI;
    double dv = 2.0 * M_PI / 80;
    double du = M_PI / 32;

    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = 180.0;
    pointLight.m[2][1] = 180.0;
    pointLight.m[3][1] = 180.0;
    pointLight.m[4][1] = 1.0;

    //points to draw sphere. the buffer here has same idea as in Torus'
    dmatrix_t P[3];
    for (int i = 0; i < 3; i++)
    {
        dmat_alloc(&P[i], 4, 1);
    }
    dmatrix_t x, y, z, buffer;
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    buffer.m[1][1] = radiumphere * cos(v) * sin(u);
    buffer.m[2][1] = radiumphere * sin(v) * sin(u);
    buffer.m[3][1] = radiumphere * cos(u);
    buffer.m[4][1] = 1;
    buffer = *perspective_projection(dmat_mult(&C, &buffer));


    /* Notice: the sphere's u, that is every slices' euqation, only need run a half */
    double delta = 2;
    double start = 0;
    for (v = (start) * M_PI; v+dv <= (delta+start) * M_PI; v += dv)
    {
        for (u = 0.0; u+du <= M_PI; u += du)
        {
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
            z.m[3][1] = radiumphere * sin(u);
            z.m[4][1] = 1;
            buffer.m[1][1] = radiumphere * cos(v + dv) * sin(u + du);
            buffer.m[2][1] = radiumphere * sin(v + dv) * sin(u + du);
            buffer.m[3][1] = radiumphere * cos(u + du);
            buffer.m[4][1] = 1;

            //2.light and calculation of brightness;
            //calcuate vectors based on surface x-y-z-buffer
            dmatrix_t vecYX;dmat_alloc(&vecYX, 4, 1);vecYX = *dmat_sub(&x, &y);//x-y
            dmatrix_t vecXZ;dmat_alloc(&vecXZ, 4, 1);vecXZ = *dmat_sub(&z, &x);
            dmatrix_t vecYZ;dmat_alloc(&vecYZ, 4, 1);vecYZ = *dmat_sub(&z, &y);
            dmatrix_t vecZbuffer;dmat_alloc(&vecZbuffer, 4, 1);vecZbuffer = *dmat_sub(&buffer, &z);

            //calcuate two vectors: xyz, yzbuffer
            dmatrix_t vecNorm_xyz; dmat_alloc(&vecNorm_xyz, 4, 1);vecNorm_xyz = *crossProduct(vecXZ, vecYX);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            dmatrix_t vecNorm_yzb; dmat_alloc(&vecNorm_yzb, 4, 1);vecNorm_yzb = *crossProduct(vecZbuffer, vecYZ);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);

            //calcuate center of xyz, yzb
            dmatrix_t pointCenter_xyz;dmat_alloc(&pointCenter_xyz,4,1);pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1])/3;pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1])/3;pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1])/3;pointCenter_xyz.m[4][1] = 1;

            dmatrix_t pointCenter_yzb;dmat_alloc(&pointCenter_yzb,4,1);pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1])/3;pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1])/3;pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1])/3;pointCenter_yzb.m[4][1] = 1;

            //calculate two light vector
            dmatrix_t vecLight_xyz;dmat_alloc(&vecLight_xyz, 4, 1);vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            dmatrix_t vecLight_yzb;dmat_alloc(&vecLight_yzb, 4, 1);vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);

            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);

            //calculate cos of two vector X light
            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = -brightness_xyz;
            //double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
            double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = -brightness_yzb;


            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;


            //3.draw centre to light part and normal vector(optional)
            //xyz
            dmatrix_t pointNormEnd_xyz; dmat_alloc(&pointNormEnd_xyz, 4, 1);pointNormEnd_xyz = *dmat_add(&pointCenter_xyz, &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&C, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&C, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);

            //yzb
            dmatrix_t pointNormEnd_yzb; dmat_alloc(&pointNormEnd_yzb, 4, 1);pointNormEnd_yzb = *dmat_add(&pointCenter_yzb, &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&C, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&C, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&C, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1]);


            //4. draw the meshes.
            x = *perspective_projection(dmat_mult(&C, &x));
            y = *perspective_projection(dmat_mult(&C, &y));
            z = *perspective_projection(dmat_mult(&C, &z));
            buffer = *perspective_projection(dmat_mult(&C, &buffer));

            //Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
            //Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
            //Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
            //5.fill with color
            P[0] = x;P[1] = y;P[2] = z;
            XFillConvexPolygon(d, w, s, P, 3, r * brightness_xyz, g * brightness_xyz, b * brightness_xyz);
            P[0] = buffer;
            XFillConvexPolygon(d, w, s, P, 3, r * brightness_yzb, g * brightness_yzb, b * brightness_yzb);
        }
    }
    //7.clean
    free_dmatrix(x.m,1,x.l,1,x.c) ;
    free_dmatrix(y.m,1,y.l,1,y.c) ;
    free_dmatrix(z.m,1,z.l,1,z.c) ;
    free_dmatrix(buffer.m,1,buffer.l,1,buffer.c) ;
    free_dmatrix(pointLight.m,1,pointLight.l,1,pointLight.c) ;
}
void cone(unsigned int red, unsigned int g, unsigned int b)
{
    //matrix designed for cone;
    dmatrix_t D = C;

    //D = *rotation('z', M_PI/12, D);
    //D = *rotation('x', 12*M_PI/12, D);
    //D = *rotation('y', -1*M_PI/12, D);
    D = *translation(0, 0, -480, D);


    //information of sphere
    double height = 400;
    double r = 30;

    //range of sphere
    double theta = 0.0 * M_PI;
    double u = 0.0;
    double dtheta = M_PI / 64;
    double du = height / 400;


    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = 180.0;
    pointLight.m[2][1] = 180.0;
    pointLight.m[3][1] = 180.0;
    pointLight.m[4][1] = 1.0;

    //points to draw sphere. the buffer here has same idea as in Torus'
    dmatrix_t P[3];
    for (int i = 0; i < 3; i++)
    {
        dmat_alloc(&P[i], 4, 1);
    }
    dmatrix_t x, y, z, buffer;
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    buffer.m[1][1] = r*cos(theta)*(height-u)/height;
    buffer.m[2][1] = r*sin(theta)*(height-u)/height;
    buffer.m[3][1] = u;
    buffer.m[4][1] = 1;
    buffer = *perspective_projection(dmat_mult(&D, &buffer));


    /* Notice: the sphere's u, that is every slices' euqation, only need run a half */
    double delta = 2;
    double start = 0;
    for  (u = height; u+du >= 0; u -= du)
    {
        for(theta = (start) * M_PI; theta <= (delta+start) * M_PI; theta += dtheta)
        {
            x.m[1][1] = r*cos(theta)*(height-u)/height;
            x.m[2][1] = r*sin(theta)*(height-u)/height;
            x.m[3][1] = u;
            x.m[4][1] = 1;
            y.m[1][1] = r*cos(theta+dtheta)*(height-u)/height;
            y.m[2][1] = r*sin(theta+dtheta)*(height-u)/height;
            y.m[3][1] = u;
            y.m[4][1] = 1;
            z.m[1][1] = r*cos(theta)*(height-u-du)/height;
            z.m[2][1] = r*sin(theta)*(height-u-du)/height;
            z.m[3][1] = u-du;
            z.m[4][1] = 1;
            buffer.m[1][1] = r*cos(theta+dtheta)*(height-u-du)/height;
            buffer.m[2][1] = r*sin(theta+dtheta)*(height-u-du)/height;
            buffer.m[3][1] = u-du;
            buffer.m[4][1] = 1;

            //2.light and calculation of brightness;
            //calcuate vectors based on surface x-y-z-buffer
            dmatrix_t vecYX;dmat_alloc(&vecYX, 4, 1);vecYX = *dmat_sub(&x, &y);//x-y
            dmatrix_t vecXZ;dmat_alloc(&vecXZ, 4, 1);vecXZ = *dmat_sub(&z, &x);
            dmatrix_t vecYZ;dmat_alloc(&vecYZ, 4, 1);vecYZ = *dmat_sub(&z, &y);
            dmatrix_t vecZbuffer;dmat_alloc(&vecZbuffer, 4, 1);vecZbuffer = *dmat_sub(&buffer, &z);

            //calcuate two vectors: xyz, yzbuffer
            dmatrix_t vecNorm_xyz; dmat_alloc(&vecNorm_xyz, 4, 1);vecNorm_xyz = *crossProduct(vecXZ, vecYX);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            dmatrix_t vecNorm_yzb; dmat_alloc(&vecNorm_yzb, 4, 1);vecNorm_yzb = *crossProduct(vecZbuffer, vecYZ);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);

            //calcuate center of xyz, yzb
            dmatrix_t pointCenter_xyz;dmat_alloc(&pointCenter_xyz,4,1);pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1])/3;pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1])/3;pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1])/3;pointCenter_xyz.m[4][1] = 1;

            dmatrix_t pointCenter_yzb;dmat_alloc(&pointCenter_yzb,4,1);pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1])/3;pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1])/3;pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1])/3;pointCenter_yzb.m[4][1] = 1;

            //calculate two light vector
            dmatrix_t vecLight_xyz;dmat_alloc(&vecLight_xyz, 4, 1);vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            dmatrix_t vecLight_yzb;dmat_alloc(&vecLight_yzb, 4, 1);vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);

            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);

            //calculate cos of two vector X light
            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = -brightness_xyz;
            //double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
            double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = -brightness_yzb;


            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;


            //3.draw centre to light part and normal vector(optional)
            //xyz
            dmatrix_t pointNormEnd_xyz; dmat_alloc(&pointNormEnd_xyz, 4, 1);pointNormEnd_xyz = *dmat_add(&pointCenter_xyz, &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&D, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&D, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&D, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);

            //yzb
            dmatrix_t pointNormEnd_yzb; dmat_alloc(&pointNormEnd_yzb, 4, 1);pointNormEnd_yzb = *dmat_add(&pointCenter_yzb, &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&D, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&D, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&D, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1]);


            //4. draw the meshes.
            x = *perspective_projection(dmat_mult(&D, &x));
            y = *perspective_projection(dmat_mult(&D, &y));
            z = *perspective_projection(dmat_mult(&D, &z));
            buffer = *perspective_projection(dmat_mult(&D, &buffer));

            //Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
            //Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
            //Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
            //5.fill with color
            P[0] = x;P[1] = y;P[2] = z;
            XFillConvexPolygon(d, w, s, P, 3, red * brightness_xyz, g * brightness_xyz, b * brightness_xyz);
            P[0] = buffer;
            XFillConvexPolygon(d, w, s, P, 3, red * brightness_yzb, g * brightness_yzb, b * brightness_yzb);
        }
    }
    //7.clean
    free_dmatrix(x.m,1,x.l,1,x.c) ;
    free_dmatrix(y.m,1,y.l,1,y.c) ;
    free_dmatrix(z.m,1,z.l,1,z.c) ;
    free_dmatrix(buffer.m,1,buffer.l,1,buffer.c) ;
    free_dmatrix(pointLight.m,1,pointLight.l,1,pointLight.c) ;
}
void reversedcone(unsigned int red, unsigned int g, unsigned int b)
{
    //matrix designed for cone;
    dmatrix_t D = C;

    //D = *rotation('z', M_PI/12, D);
    //D = *rotation('x', 12*M_PI/12, D);
    D = *rotation('y', 12*M_PI/12, D);
    D = *translation(0, 0, -480, D);


    //information of sphere
    double height = 400;
    double r = 30;

    //range of sphere
    double theta = 0.0 * M_PI;
    double u = 0.0;
    double dtheta = M_PI / 64;
    double du = height / 400;


    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = 180.0;
    pointLight.m[2][1] = 180.0;
    pointLight.m[3][1] = 180.0;
    pointLight.m[4][1] = 1.0;

    //points to draw sphere. the buffer here has same idea as in Torus'
    dmatrix_t P[3];
    for (int i = 0; i < 3; i++)
    {
        dmat_alloc(&P[i], 4, 1);
    }
    dmatrix_t x, y, z, buffer;
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&buffer, 4, 1);
    buffer.m[1][1] = r*cos(theta)*(height-u)/height;
    buffer.m[2][1] = r*sin(theta)*(height-u)/height;
    buffer.m[3][1] = u;
    buffer.m[4][1] = 1;
    buffer = *perspective_projection(dmat_mult(&D, &buffer));


    /* Notice: the sphere's u, that is every slices' euqation, only need run a half */
    double delta = 2;
    double start = 0;
    for  (u = height; u+du >= 0; u -= du)
    {
        for(theta = (start) * M_PI; theta <= (delta+start) * M_PI; theta += dtheta)
        {
            x.m[1][1] = r*cos(theta)*(height-u)/height;
            x.m[2][1] = r*sin(theta)*(height-u)/height;
            x.m[3][1] = u;
            x.m[4][1] = 1;
            y.m[1][1] = r*cos(theta+dtheta)*(height-u)/height;
            y.m[2][1] = r*sin(theta+dtheta)*(height-u)/height;
            y.m[3][1] = u;
            y.m[4][1] = 1;
            z.m[1][1] = r*cos(theta)*(height-u-du)/height;
            z.m[2][1] = r*sin(theta)*(height-u-du)/height;
            z.m[3][1] = u-du;
            z.m[4][1] = 1;
            buffer.m[1][1] = r*cos(theta+dtheta)*(height-u-du)/height;
            buffer.m[2][1] = r*sin(theta+dtheta)*(height-u-du)/height;
            buffer.m[3][1] = u-du;
            buffer.m[4][1] = 1;

            //2.light and calculation of brightness;
            //calcuate vectors based on surface x-y-z-buffer
            dmatrix_t vecYX;dmat_alloc(&vecYX, 4, 1);vecYX = *dmat_sub(&x, &y);//x-y
            dmatrix_t vecXZ;dmat_alloc(&vecXZ, 4, 1);vecXZ = *dmat_sub(&z, &x);
            dmatrix_t vecYZ;dmat_alloc(&vecYZ, 4, 1);vecYZ = *dmat_sub(&z, &y);
            dmatrix_t vecZbuffer;dmat_alloc(&vecZbuffer, 4, 1);vecZbuffer = *dmat_sub(&buffer, &z);

            //calcuate two vectors: xyz, yzbuffer
            dmatrix_t vecNorm_xyz; dmat_alloc(&vecNorm_xyz, 4, 1);vecNorm_xyz = *crossProduct(vecXZ, vecYX);
            vecNorm_xyz = *dmat_normalize(&vecNorm_xyz);vecNorm_xyz = *dmat_scalar_mult(&vecNorm_xyz, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            dmatrix_t vecNorm_yzb; dmat_alloc(&vecNorm_yzb, 4, 1);vecNorm_yzb = *crossProduct(vecZbuffer, vecYZ);
            vecNorm_yzb = *dmat_normalize(&vecNorm_yzb);vecNorm_yzb = *dmat_scalar_mult(&vecNorm_yzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);

            //calcuate center of xyz, yzb
            dmatrix_t pointCenter_xyz;dmat_alloc(&pointCenter_xyz,4,1);pointCenter_xyz.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1])/3;pointCenter_xyz.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1])/3;pointCenter_xyz.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1])/3;pointCenter_xyz.m[4][1] = 1;

            dmatrix_t pointCenter_yzb;dmat_alloc(&pointCenter_yzb,4,1);pointCenter_yzb.m[1][1] = (z.m[1][1] + y.m[1][1] + buffer.m[1][1])/3;pointCenter_yzb.m[2][1] = (z.m[2][1] + y.m[2][1] + buffer.m[2][1])/3;pointCenter_yzb.m[3][1] = (z.m[3][1] + y.m[3][1] + buffer.m[3][1])/3;pointCenter_yzb.m[4][1] = 1;

            //calculate two light vector
            dmatrix_t vecLight_xyz;dmat_alloc(&vecLight_xyz, 4, 1);vecLight_xyz = *dmat_sub(&pointLight, &pointCenter_xyz);
            dmatrix_t vecLight_yzb;dmat_alloc(&vecLight_yzb, 4, 1);vecLight_yzb = *dmat_sub(&pointLight, &pointCenter_yzb);

            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);

            //calculate cos of two vector X light
            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = -brightness_xyz;
            //double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = 0;
            double brightness_yzb = angle(vecNorm_yzb,vecLight_yzb);if (brightness_yzb < 0)brightness_yzb = -brightness_yzb;


            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;


            //3.draw centre to light part and normal vector(optional)
            //xyz
            dmatrix_t pointNormEnd_xyz; dmat_alloc(&pointNormEnd_xyz, 4, 1);pointNormEnd_xyz = *dmat_add(&pointCenter_xyz, &vecNorm_xyz);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&D, &pointNormEnd_xyz));
            pointCenter_xyz = *perspective_projection(dmat_mult(&D, &pointCenter_xyz));
            vecLight_xyz = *perspective_projection(dmat_mult(&D, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);

            //yzb
            dmatrix_t pointNormEnd_yzb; dmat_alloc(&pointNormEnd_yzb, 4, 1);pointNormEnd_yzb = *dmat_add(&pointCenter_yzb, &vecNorm_yzb);//pointNormEnd_yzb = pointCenter_xyz+Normalvector
            pointNormEnd_yzb = *perspective_projection(dmat_mult(&D, &pointNormEnd_yzb));
            pointCenter_yzb = *perspective_projection(dmat_mult(&D, &pointCenter_yzb));
            vecLight_yzb = *perspective_projection(dmat_mult(&D, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_yzb.m[1][1], pointCenter_yzb.m[2][1], pointNormEnd_yzb.m[1][1], pointNormEnd_yzb.m[2][1]);


            //4. draw the meshes.
            x = *perspective_projection(dmat_mult(&D, &x));
            y = *perspective_projection(dmat_mult(&D, &y));
            z = *perspective_projection(dmat_mult(&D, &z));
            buffer = *perspective_projection(dmat_mult(&D, &buffer));

            //Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
            //Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
            //Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
            //5.fill with color
            P[0] = x;P[1] = y;P[2] = z;
            XFillConvexPolygon(d, w, s, P, 3, red * brightness_xyz, g * brightness_xyz, b * brightness_xyz);
            P[0] = buffer;
            XFillConvexPolygon(d, w, s, P, 3, red * brightness_yzb, g * brightness_yzb, b * brightness_yzb);
        }
    }
    //7.clean
    free_dmatrix(x.m,1,x.l,1,x.c) ;
    free_dmatrix(y.m,1,y.l,1,y.c) ;
    free_dmatrix(z.m,1,z.l,1,z.c) ;
    free_dmatrix(buffer.m,1,buffer.l,1,buffer.c) ;
    free_dmatrix(pointLight.m,1,pointLight.l,1,pointLight.c) ;
}

//cone, using square to render.
void cone2(unsigned int red, unsigned int g, unsigned int b)
{
    //matrix designed for cone;
    dmatrix_t D = C;

    //D = *scalarization(1.5, 1.5, 1.5,D);

    //D = *rotation('z', M_PI/12, D);
    //D = *rotation('x', 12*M_PI/12, D);
    //D = *rotation('y', M_PI/12, D);
    //D = *translation(0, 0, -180, D);


    //information of sphere
    double height = 100;
    double r = 50;

    //range of sphere
    double theta = 0.0 * M_PI;
    double u = 0.0;
    double dtheta = M_PI / 8;
    double du = height / 50;


    //light
    dmatrix_t pointLight;
    dmat_alloc(&pointLight, 4, 1);
    pointLight.m[1][1] = 180.0;
    pointLight.m[2][1] = 180.0;
    pointLight.m[3][1] = 180.0;
    pointLight.m[4][1] = 1.0;

    //points to draw sphere. the buffer here has same idea as in Torus'
    dmatrix_t P[4];
    for (int i = 0; i < 4; i++)
    {
        dmat_alloc(&P[i], 4, 1);
    }
    dmatrix_t x, y, z, buffer;
    dmat_alloc(&x, 4, 1);
    dmat_alloc(&y, 4, 1);
    dmat_alloc(&z, 4, 1);
    dmat_alloc(&buffer, 4, 1);

    /* Notice: the sphere's u, that is every slices' euqation, only need run a half */
    double delta = 2;
    double start = 0;
    for  (u = height; u-du >= 0; u -= du)
    {
        for(theta = (start) * M_PI; theta <= (delta+start) * M_PI; theta += dtheta)
        {
            x.m[1][1] = r*cos(theta)*(height-u)/height;
            x.m[2][1] = r*sin(theta)*(height-u)/height;
            x.m[3][1] = u;
            x.m[4][1] = 1;
            y.m[1][1] = r*cos(theta+dtheta)*(height-u)/height;
            y.m[2][1] = r*sin(theta+dtheta)*(height-u)/height;
            y.m[3][1] = u;
            y.m[4][1] = 1;
            z.m[1][1] = r*cos(theta)*(height-u-du)/height;
            z.m[2][1] = r*sin(theta)*(height-u-du)/height;
            z.m[3][1] = u-du;
            z.m[4][1] = 1;
            buffer.m[1][1] = r*cos(theta+dtheta)*(height-u-du)/height;
            buffer.m[2][1] = r*sin(theta+dtheta)*(height-u-du)/height;
            buffer.m[3][1] = u-du;
            buffer.m[4][1] = 1;

            //2.light and calculation of brightness;
            //calcuate vectors based on surface x-y-z-buffer
            dmatrix_t vecYX;dmat_alloc(&vecYX, 4, 1);vecYX = *dmat_sub(&x, &y);//x-y
            dmatrix_t vecXZ;dmat_alloc(&vecXZ, 4, 1);vecXZ = *dmat_sub(&z, &x);
            dmatrix_t vecYZ;dmat_alloc(&vecYZ, 4, 1);vecYZ = *dmat_sub(&z, &y);
            dmatrix_t vecZbuffer;dmat_alloc(&vecZbuffer, 4, 1);vecZbuffer = *dmat_sub(&buffer, &z);

            //calcuate two vectors: xyz, yzbuffer
            dmatrix_t vecNorm_xyzb; dmat_alloc(&vecNorm_xyzb, 4, 1);vecNorm_xyzb = *crossProduct(vecXZ, vecYX);
            vecNorm_xyzb = *dmat_normalize(&vecNorm_xyzb);vecNorm_xyzb = *dmat_scalar_mult(&vecNorm_xyzb, 100); //normailize it and times 10 could make it looks clear when visualize normalvector.

            //std::cout<<"yx:"<<std::endl;write_dmatrix(&vecYX);
            //std::cout<<"xz:"<<std::endl;write_dmatrix(&vecXZ);

            //calcuate center of xyz, yzb
            dmatrix_t pointCenter_xyzb;dmat_alloc(&pointCenter_xyzb,4,1);pointCenter_xyzb.m[1][1] = (x.m[1][1] + y.m[1][1] + z.m[1][1] + buffer.m[1][1])/4;pointCenter_xyzb.m[2][1] = (x.m[2][1] + y.m[2][1] + z.m[2][1] + buffer.m[2][1])/4;pointCenter_xyzb.m[3][1] = (x.m[3][1] + y.m[3][1] + z.m[3][1] + buffer.m[3][1])/4;pointCenter_xyzb.m[4][1] = 1;

            //calculate two light vector
            dmatrix_t vecLight_xyzb;dmat_alloc(&vecLight_xyzb, 4, 1);vecLight_xyzb = *dmat_sub(&pointLight, &vecLight_xyzb);

            //std::cout<<"\npointCenter_xyz:"<<std::endl; write_dmatrix(&pointCenter_xyz);
            //std::cout<<"\nsourceLight:"<<std::endl;write_dmatrix(&pointLight);
            //std::cout<<"\n***vecLight_xyz:"<<std::endl;write_dmatrix(&vecLight_xyz);
            //std::cout<<"\n***vecNorm_xyz:"<<std::endl; write_dmatrix(&vecNorm_xyz);

            //calculate cos of two vector X light
            //double brightness_xyz = angle(vecNorm_xyz,vecLight_xyz);if (brightness_xyz < 0)brightness_xyz = 0;
            double brightness = angle(vecNorm_xyzb,vecLight_xyzb);if (brightness < 0)brightness = -brightness;

            //std::cout<<"\nbrightness: "<<brightness<<"\n======="<<std::endl;


            //3.draw centre to light part and normal vector(optional)
            //xyz
            dmatrix_t pointNormEnd_xyz; dmat_alloc(&pointNormEnd_xyz, 4, 1);pointNormEnd_xyz = *dmat_add(&pointCenter_xyzb, &vecNorm_xyzb);//pointNormEnd_xyz = pointCenter_xyz+Normalvector
            pointNormEnd_xyz = *perspective_projection(dmat_mult(&D, &pointNormEnd_xyz));
            pointCenter_xyzb = *perspective_projection(dmat_mult(&D, &pointCenter_xyzb));
            vecLight_xyzb = *perspective_projection(dmat_mult(&D, &pointLight));
            //light to point, and normal vector
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], vecLight_xyz.m[1][1], vecLight_xyz.m[2][1]);
            //Line(pointCenter_xyz.m[1][1], pointCenter_xyz.m[2][1], pointNormEnd_xyz.m[1][1], pointNormEnd_xyz.m[2][1]);

            //4. draw the meshes.
            x = *perspective_projection(dmat_mult(&D, &x));
            y = *perspective_projection(dmat_mult(&D, &y));
            z = *perspective_projection(dmat_mult(&D, &z));
            buffer = *perspective_projection(dmat_mult(&D, &buffer));
//            Line(x.m[1][1],x.m[2][1],y.m[1][1],y.m[2][1]);
//            Line(x.m[1][1],x.m[2][1],z.m[1][1],z.m[2][1]);
//            Line(y.m[1][1],y.m[2][1],z.m[1][1],z.m[2][1]);
            //5.fill with color
            P[0] = x;P[1] = y;P[2] = z;  P[3] = buffer;
            XFillConvexPolygon(d, w, s, P, 4, red * brightness, g * brightness, b * brightness);
        }
    }
    //7.clean
    free_dmatrix(x.m,1,x.l,1,x.c) ;
    free_dmatrix(y.m,1,y.l,1,y.c) ;
    free_dmatrix(z.m,1,z.l,1,z.c) ;
    free_dmatrix(buffer.m,1,buffer.l,1,buffer.c) ;
    free_dmatrix(pointLight.m,1,pointLight.l,1,pointLight.c) ;
}


void Draw()
{
    hiddenremover(); //apply hidden algorithm
    cone(255, 248, 220);
    torus(220, 220, 220);
    sphere(237, 189, 101);
    reversedcone(255, 248, 220);
}

void OnDisplay()
{
    memset(frame, 255, windowW * windowH * 3);
    //main print function
    Draw();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(windowW, windowH, GL_RGB, GL_UNSIGNED_BYTE, (GLubyte *) frame);
    glFlush();
    glutSwapBuffers();
}


