#include <iostream>
#include "draw.h"

void idle()
{
    rotation('x', 0.5);
    rotation('y', 0.5);
    rotation('z', 0.5);
    angle_total += 0.5;
    //scalarization(0.01, 0.01, 0.01);
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

    scalarization(1, 1, 1);
    rotation('x', (2 * M_PI / 12));
    rotation('y', (1 * M_PI / 12));
    rotation('z', (7 * M_PI / 12));

    //-- run the program
    glutDisplayFunc(OnDisplay);
    glutKeyboardFunc(OnKeyboard);
    //glutIdleFunc(idle);
    glutMainLoop();
    delete_dmatrix(&C);
    return 0;
}
