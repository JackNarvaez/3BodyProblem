#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

void Acceleration(const double q0[][3], const double mass[], double a[][3]);

void PEFRL(const double q0[][6], const double mass[], const int &n, const double &dt, double q[][3][6]);

void rel_vec(const double Center_Body[], const double Outer_Body[], double Relative_Position[], const int &n);

double Energy(double q[][6], double mass[]);

#endif