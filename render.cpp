#include <iostream>
#include "draw.h"

const double rate = 0.02;

void idle()
{
    rotation('x', +rate);
    rotation('y', -rate);
    rotation('z', +rate);
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
    rotation('x', (10  * M_PI) / 24);
    rotation('y', (21 * M_PI) / 24);
    rotation('z', (13 * M_PI) / 24);

    glutDisplayFunc(OnDisplay);
    glutKeyboardFunc(OnKeyboard);
    glutIdleFunc(idle);
    glutMainLoop();
    delete_dmatrix(&C);
    return 0;
}
