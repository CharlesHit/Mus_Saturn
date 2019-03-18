#pragma once

const int windowW = 1080;
const int windowH = 1080;
const double ASPECT = (double)windowW / (double)windowH;

const double THETA = 90.0;

//the camera's eye and gaze points.
const double Ex = 150;
const double Ey = 150;
const double Ez = -150;

const double Gx = 0.0;
const double Gy = 0.0;
const double Gz = 0.0;

const double UPx = 0.0;
const double UPy = 0.0;
const double UPz = 1;

//const double direction[3] = { UPx,UPy,UPz }; //the direction vector for camera

const double FP = 50; //the far point for camera
const double NP = 5; //the near point for camera
