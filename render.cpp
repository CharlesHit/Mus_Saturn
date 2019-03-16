#include <iostream>
#include "draw.h"

void idle()
{
    rotation('x', 0.001);
    rotation('y', 0.001);
    rotation('z', 0.001);
    angle_total += 0.001;
    //scalarization(0.01, 0.01, 0.01);
    OnDisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(windowW, windowH);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Mu He");

    glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel(GL_FLAT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    cameraInitialization();

    scalarization(0.8, 0.8, 0.8);

    rotation('y', 1.505000+3 * M_PI / 12);
    rotation('x', 1.505000+5 * M_PI / 12);
    rotation('z', 1.505000);

    //-- run the program
    glutDisplayFunc(OnDisplay);
    glutKeyboardFunc(OnKeyboard);
    glutIdleFunc(idle);
    glutMainLoop();
    delete_dmatrix(&C);
    return 0;
}
