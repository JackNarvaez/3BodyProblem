#include <iostream>
#include <cmath>
#include "Integrator.h"

double Energy(double q[][6], double mass){
    /*---------------------------------------------------------------------------------------------
    Energy:
    Calculate the total energy of 3 particles interacting gravitationally.
    -----------------------------------------------------------------------------------------------
    Arguments:
        q   :   Array with the position and velocity of the particles with the format
                q[i] = [xi, yi, zi, vxi, vyi, vzi]
        mass:   Array with the masses of the particles.
                m = [m1, m2, m3]
    -----------------------------------------------------------------------------------------------
    Returns:
    E = Total energy of the system
    ---------------------------------------------------------------------------------------------*/
    double speed2[3];
    for (int ii=0; ii<3; ii++){
        speed2[ii] = pow(q[ii][3],2)+pow(q[ii][4],2)+pow(q[ii][5],2);
    }
    double rel_position[3][3], r[3];
    rel_vec(q[1],q[0],rel_position[0]);
    rel_vec(q[2],q[0],rel_position[0]);
    rel_vec(q[2],q[1],rel_position[0]);

    for (int ii=0, ii < 3; ii++){
        r[ii] = sqrt(pow(rel_position[ii][0],2)+pow(rel_position[ii][1],2)+pow(rel_position[ii][2],2));
    }
    double U[3] = {mass[0]*mass[1]/r[0], mass[0]*mass[2]/r[1], mass[1]*mass[2]/r[2]};
    double E=0;
    for (int ii=0; ii < 3; ii++){
        E+=0.5*speed2[ii]*mass[ii]-2*G*U[ii];
    }    
    return E;
}

double Ang_momentum(double q[][6], double mass[]){
    /*---------------------------------------------------------------------------------------------
    Constants:
    Calculates the total angular momentum of 3 particles interacting gravitationally.
    -----------------------------------------------------------------------------------------------
    Arguments:
        q   :   Array with the position and velocity of the particles with the format
                q[i] = [xi, yi, zi, vxi, vyi, vzi]
        mass:   Array with the masses of the particles.
                m = [m1, m2, m3]
    -----------------------------------------------------------------------------------------------
    Returns:
        L = Magnitude of the total angular momentum the system
    ---------------------------------------------------------------------------------------------*/
    double ang_mom_vector[3] = {0,0,0};
    for (int ii=0; ii<3; ii++){
        ang_mom_vector[0] += mass[ii]*(q[ii][1]*q[ii][5]-q[ii][2]*q[ii][4]);
        ang_mom_vector[1] += mass[ii]*(q[ii][2]*q[ii][3]-q[ii][0]*q[ii][5]);
        ang_mom_vector[2] += mass[ii]*(q[ii][0]*q[ii][4]-q[ii][1]*q[ii][3]);
    }    
    return sqrt(pow(ang_mom_vector[0],2)+pow(ang_mom_vector[1],2)+pow(ang_mom_vector[2],2));
}