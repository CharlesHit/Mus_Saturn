#pragma once

const int windowW = 1000;
const int windowH = 800;
const double ASPECT = (double)windowW / (double)windowH;

const double THETA = 90.0;

//the camera's eye and gaze points.
const double Ex = 180;
const double Ey = 180;
const double Ez = 180;

const double E[3] = { Ex, Ey, Ez };

const double Gx = 0.0;
const double Gy = 0.0;
const double Gz = 0.0;

const double G[3] = { Gx, Gy, Gz };

const double UPx = 0.0;
const double UPy = 1.0;
const double UPz = 0.0;

//const double direction[3] = { UPx,UPy,UPz }; //the direction vector for camera

const double FP = 50; //the far point for camera
const double NP = 5; //the near point for camera
