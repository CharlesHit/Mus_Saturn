#include <iostream>
#include "draw.h"
#include "math.h"

const double delta_theta = M_PI/32, rad = 30;
double theta = 0;

void idle()
{
    const double rate = 0.02;
    theta = theta + delta_theta;

    rotation('x', +rate);
    rotation('y', +rate);
    rotation('z', +rate);
    translation(rad * cos(theta), rad * sin(theta), 0);
    angle_total += rate;
    OnDisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(windowW, windowH);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Mu He - Mu's Saturn");

    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_FLAT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    cameraInitialization();

//    scalarization(1, 1, 1);
    scalarization(0.2, 0.2, 0.2);
//    translation(10, 10, 10);
//    rotation('x', (10  * M_PI) / 24);
//    rotation('y', (21 * M_PI) / 24);
//    rotation('z', (13 * M_PI) / 24);


    glutDisplayFunc(OnDisplay);
    glutKeyboardFunc(OnKeyboard);
    glutIdleFunc(idle);
    glutMainLoop();
    delete_dmatrix(&C);

    return 0;
}
