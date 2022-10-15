// Simulation of a system of three bodies interacting with each other gravitationally.
// Evolution is performed by PEFRL integrator.

#include <cmath>
#include "Integrator.h"

const double G{4*pow(M_PI,2)}; // Universal gravitational constant.

void Acceleration(const double q0[][3], const double mass[], double a[][3]){
    /*---------------------------------------------------------------------------------------------
    Calculates the gravitational acceleration for each body in the 3-body system.
    -----------------------------------------------------------------------------------------------
    Arguments:
    q0: Array with the initial condition data:
        q0[0] = particle 1
        q0[1] = particle 2
        q0[2] = particle 3, where:
        q0[i] = [x_i, y_i, z_i] for i=0,1,2
    mass: masses of the particles
        mass = [m1, m2, m3]
    a:  Array that stores the gravitational acceleration.
    -----------------------------------------------------------------------------------------------
    Fills the values of the acceleration vector, as follows:
        a[i] = [ax_i, ay_i, az_i] for i=0,1,2
    ---------------------------------------------------------------------------------------------*/
    double Deltaxyz[9]{0}, r[3]{0};
    for (int jj = 0; jj < 3; jj++){
        for (int ii = 1; ii < 3; ii++)
            Deltaxyz[ii*3+jj-3]=q0[0][jj]-q0[ii][jj];
        Deltaxyz[6+jj] = q0[1][jj]-q0[2][jj];
    }
    for (int ii=0; ii < 3; ii++){
        r[ii]= sqrt(Deltaxyz[ii*3]*Deltaxyz[ii*3] + Deltaxyz[ii*3+1]*Deltaxyz[ii*3+1] 
                    +Deltaxyz[ii*3+2]*Deltaxyz[ii*3+2]);
    }
    for (int jj = 0; jj < 3; jj++){
        a[0][jj] = -G*(Deltaxyz[jj]*mass[1]/(r[0]*r[0]*r[0]) + Deltaxyz[3+jj] * mass[2]/(r[1]*r[1]*r[1]));
        a[1][jj] = G*(Deltaxyz[jj]*mass[0]/(r[0]*r[0]*r[0]) - Deltaxyz[6+jj] * mass[2]/(r[2]*r[2]*r[2]));
        a[2][jj] = G*(Deltaxyz[3+jj]*mass[0]/(r[1]*r[1]*r[1]) + Deltaxyz[6+jj] * mass[1]/(r[2]*r[2]*r[2]));
    }
}

void PEFRL(const double q0[][6], const double mass[], const int &n, const double &dt, double q[][3][6]){
    /*---------------------------------------------------------------------------------------------
    Position Extended Forest-Ruth Like (PEFRL) method for time evolution.
    -----------------------------------------------------------------------------------------------
    Arguments:
    q0  :   Values of ODEs system at t=t0. The format is:
            q[i]:   [x_i, y_i, z_i, vx_i, vy_i, vz_i] for i = 0,1,2.
    mass:   Masses of the particles. The format is:
            mass = [m1, m2, m3]
    n   :   Number of steps in the grid.
    dt  :   Stepsize for the iteration.
    -----------------------------------------------------------------------------------------------
    Fills q array with the components of the solution since t0 to tf=ndt.
    ---------------------------------------------------------------------------------------------*/
    // Arrays that store acceleration, position, and velocity, respectively, in intermediate steps.
    double F[3][3]{0}, X[3][3]{0}, V[3][3]{0};
    // Initial condition
    for (int ii = 0; ii < 3; ii++){
        for (int jj = 0; jj < 6; jj++){
            q[0][ii][jj] = q0[ii][jj];
        }
    }
    // PEFRL Parameters
    double xi{0.1786178958448091e+0};
    double gamma{-0.2123418310626054e+0};
    double chi{-0.6626458266981849e-1};
    // Main loop
    for (int ll = 0; ll < n-1; ll++){
        for (int ii = 0; ii < 3; ii++){
            for (int jj = 0; jj < 3; jj++){
                X[ii][jj] = q[ll][ii][jj] + xi*dt*q[ll][ii][3+jj];
            }
        }
        Acceleration(X, mass, F);
        for (int ii = 0; ii < 3; ii++){
            for (int jj = 0; jj < 3; jj++){
                V[ii][jj] = q[ll][ii][jj+3] + 0.5*(1.-2*gamma)*dt*F[ii][jj];
                X[ii][jj] += chi*dt*V[ii][jj];
            }
        }
        Acceleration(X, mass, F);
        for (int ii = 0; ii < 3; ii++){
            for (int jj = 0; jj < 3; jj++){
                V[ii][jj] += gamma*dt*F[ii][jj];
                X[ii][jj] += (1.-2*(chi+xi))*dt*V[ii][jj];
            }
        }
        Acceleration(X, mass, F);
        for (int ii = 0; ii < 3; ii++){
            for (int jj = 0; jj < 3; jj++){
                V[ii][jj] += gamma*dt*F[ii][jj];
                 X[ii][jj] += chi*dt*V[ii][jj];
            }
        }
        Acceleration(X, mass, F);
        for (int ii = 0; ii < 3; ii++){
            for (int jj = 0; jj < 3; jj++){
                q[ll+1][ii][jj+3] = V[ii][jj]+0.5*(1.-2*gamma)*dt*F[ii][jj];
                q[ll+1][ii][jj] = X[ii][jj]+xi*dt*q[ll+1][ii][jj+3];
            }
        }
    }
}

void rel_vec(const double Center_Body[], const double Outer_Body[], double Relative_Position[], const int &n){
    /*---------------------------------------------------------------------------------------------
    Calculates the relative vector of the outer body with respect to the center body.
    -----------------------------------------------------------------------------------------------
    Arguments:
    Center_Body :   Array with cartesian vector of Center body.
    Outer_Body  :   Array with cartesian vector of Outer body.
    n           :   Length of Arrays.
    -----------------------------------------------------------------------------------------------
    Fills Relative_Position array as follows:
        Relative_position[i] = Outer_Body[i]-Center_Body[i] for i=0,...,n
    ---------------------------------------------------------------------------------------------*/
    for (int i=0; i<n; i++){
        Relative_Position[i]=Outer_Body[i]-Center_Body[i];
    }    
}

double Energy(double q[][6], double mass[]){
    /*---------------------------------------------------------------------------------------------
    Calculates the total energy of 3 bodies interacting gravitationally.
    -----------------------------------------------------------------------------------------------
    Arguments:
    q   :   Array with the position and velocity of the particles. The format is:
                q[i] = [xi, yi, zi, vxi, vyi, vzi]
    mass:   Array with the masses of the particles.
                m = [m1, m2, m3]
    -----------------------------------------------------------------------------------------------
    Returns:
    E = Total energy of the system
    ---------------------------------------------------------------------------------------------*/
    double speed2[3]{0};
    for (int ii=0; ii<3; ii++){
        speed2[ii] = q[ii][3]*q[ii][3]+q[ii][4]*q[ii][4]+q[ii][5]*q[ii][5];
    }
    double rel_position[3][3]{0}, r[3]{0};
    rel_vec(q[1],q[0],rel_position[0], 3);
    rel_vec(q[2],q[0],rel_position[1], 3);
    rel_vec(q[2],q[1],rel_position[2], 3);
    for (int ii=0; ii < 3; ii++){
        r[ii] = sqrt(rel_position[ii][0]*rel_position[ii][0] + rel_position[ii][1]*rel_position[ii][1]
                    +rel_position[ii][2]*rel_position[ii][2]);
    }
    double U[3]{mass[0]*mass[1]/r[0], mass[0]*mass[2]/r[1], mass[1]*mass[2]/r[2]};
    double E{0};
    for (int ii=0; ii < 3; ii++){
        E+=0.5*speed2[ii]*mass[ii]-2*G*U[ii];
    }
    return E;
}