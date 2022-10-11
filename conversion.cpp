// Conversion between Orbital Elements and Cartesian State Vectors: positions and velocities.
// OEs = {a, e, i, omega, Ω, t0}
// CVs = {x, y, z, vx, vy, vz} 

#include <iostream>
#include <cmath>
#include "conversion.h"

double Halley_method(double x_0, const double &M, const double &ecc, const double &tol, const int &NMAX){
    /*---------------------------------------------------------------------------------------------
    Calculate the root, eccentric anomaly, in the Kepler's equation:
    M = E - e*sin(E),
    by Halley's method.
    -----------------------------------------------------------------------------------------------
    Arguments:
        x_0 :   Initial approximation.
        M   :   Mean anomaly.
        ecc :   Eccentricity.
        tol :   Tolerance.
        NMAX:   Maximum number of iterations.
    -----------------------------------------------------------------------------------------------
    Returns:
        E   :   Eccentric anomaly.
    ---------------------------------------------------------------------------------------------*/
    double E = std::fmod(x_0, 2*M_PI); // Find the floating-point remainder at [0, 2 pi)
    double change = 2*tol; //Difference between E_n+1 and E_n [initial value is 2*tol]
    int i = 1;  // Iteration counter
    while(std::fabs(change) > tol){
        change = 2 * (E - ecc * std::sin(E) - M) * (1. - ecc * std::cos(E))/
                    (2 * pow(1. - ecc * std::cos(E), 2) - (E - ecc * std::sin(E) - M) * ecc * std::sin(E));
        E -= change;    // Update E
        i++;
        if (i >= NMAX){
            i = NMAX+2;
            break;
        }
    }
    if (i==NMAX+2){
        std::cout <<  "Method failed after " << NMAX << "iterations." << std::endl;
        return 0;
    }
    else{
        return E;
    }
}

void op_to_coords(const double &mu, const double &a, const double &ecc, const double &i, const double &omega,
                  const double &Omega, const double &t_0, const double &t, double Cartesian_vector[]){
    /*---------------------------------------------------------------------------------------------
    Transforms from Orbital Elements to Cartesian State Vector.
    Note: Check the system of units of Arguments.
    -----------------------------------------------------------------------------------------------
    Arguments:
        mu  :   G*M, where M = central body's mass.
        a   :   Semi-major axis.
        ecc :   Eccentricity [rad].
        i   :   Inclination [rad].
        omega:  Argument of the pericenter [rad].
        Omega:  Longitude of the ascending node [rad].
        t_0 :   Epoch.
        t   :   Time.
        Cartesian_vector: Array that stores Cartesian State Vector.
    -----------------------------------------------------------------------------------------------
    Fills Cartesian_vector with the values of Cartesian State Vector, as follows:
    Cartesian_vector[] = {x, y, z, vx, vy, vz}
    ---------------------------------------------------------------------------------------------*/
    double n = sqrt(mu/(a*a*a));   // Mean motion
    double M = n * (t - t_0);   // Mean anomaly
    double E = Halley_method(0., M, ecc, 10e-10, 100);  // Eccentric anomaly
    double r = a * (1 - ecc * std::cos(E)); // Distance to the central body
    // Position vector in orbital frame
    double r_XYZ[3] = {a * (std::cos(E) - ecc), a * sqrt(1 - ecc*ecc)*std::sin(E), 0.};
    // Velocity vector in the orbital frame
    double v_XYZ[3] = {-(sqrt(mu*a)/r)*std::sin(E),(sqrt(mu*a*(1-ecc*ecc))/r)*std::cos(E),0.};
    // Transformation matrix
    double R[3][3] = {{std::cos(omega)*std::cos(Omega)-std::sin(omega)*std::sin(Omega)*std::cos(i),
                       -std::sin(omega)*std::cos(Omega)-std::cos(omega)*std::sin(Omega)*std::cos(i),
                       std::sin(Omega)*std::sin(i)},
                      {std::cos(omega)*std::sin(Omega)+std::sin(omega)*std::cos(Omega)*std::cos(i),
                       -std::sin(omega)*std::sin(Omega)+std::cos(omega)*std::cos(Omega)*std::cos(i),
                       -std::cos(Omega)*std::sin(i)},
                      {std::sin(omega)*std::sin(i),std::cos(omega)*std::sin(i),std::cos(i)}};
    // Cartesian vector in inertial frame.
    for (int l=0; l<3; l++){
        Cartesian_vector[l] = R[l][0]*r_XYZ[0]+R[l][1]*r_XYZ[1]+R[l][2]*r_XYZ[2];
        Cartesian_vector[l+3] = R[l][0]*v_XYZ[0]+R[l][1]*v_XYZ[1]+R[l][2]*v_XYZ[2];
    }
}

void coords_to_op(const double &mu, const double Cartesian_vector[], double Orbital_Elements[]){
    /*---------------------------------------------------------------------------------------------
    Transforms from Cartesian State Vector at t=0 to the Orbital Elements.
    Note: Check the system of units of Arguments.
    -----------------------------------------------------------------------------------------------
    Arguments: 
        mu      :
        Cartesian_vector:  Array with the components of the Cartesian State Vector. The format is:
                            [x, y, z, vx, vy, vz].
        Orbital_Elements:  Array that stores Orbital Elements. The format is:
                            [a, ecc, i, omega, Omega].        
    -----------------------------------------------------------------------------------------------
    Fills Orbital_Elements with the values of Orbital Elements, as follows:
        Orbital_Elements[i] = xi, where:
        x0 = a:      Semi-major axis.
        x1 = ecc:    Eccentricity [rad].
        x2 = i:      Inclination [rad].
        x3 = omega:  Argument of the pericenter [rad].
        x4 = Omega:  Longitude of the ascending node [rad].
    ---------------------------------------------------------------------------------------------*/
    // Position's norm
    double r = sqrt(Cartesian_vector[0]*Cartesian_vector[0] + Cartesian_vector[1]*Cartesian_vector[1]
                    + Cartesian_vector[2]*Cartesian_vector[2]);
    // Velocity's norm
    double v = sqrt(Cartesian_vector[3]*Cartesian_vector[3] + Cartesian_vector[4]*Cartesian_vector[4]
                    + Cartesian_vector[5]*Cartesian_vector[5]);
    // Radial velocity
    double v_r = (Cartesian_vector[0]*Cartesian_vector[3]+Cartesian_vector[1]*Cartesian_vector[4]
                    +Cartesian_vector[2]*Cartesian_vector[5])/r;
    // Angular momentum
    double h_xyz[3] = {Cartesian_vector[1]*Cartesian_vector[5]-Cartesian_vector[2]*Cartesian_vector[4],
                       Cartesian_vector[2]*Cartesian_vector[3]-Cartesian_vector[0]*Cartesian_vector[5],
                       Cartesian_vector[0]*Cartesian_vector[4]-Cartesian_vector[1]*Cartesian_vector[3]};
    // Norm of the angular momentum
    double h = sqrt(h_xyz[0]*h_xyz[0] + h_xyz[1]*h_xyz[1] + h_xyz[2]*h_xyz[2]);
    // Inclination of the orbit
    double i = std::acos(h_xyz[2]/h);
    // Line of Nodes
    double N = sqrt(h_xyz[0]*h_xyz[0] + h_xyz[1]*h_xyz[1]);
    // Longitude of ascending node
    double omega,Omega;
    if (h_xyz[0] >= 0){
        Omega = std::acos(-h_xyz[1]/N);
    }
    else{
        Omega = 2*M_PI - std::acos(-h_xyz[1]/N);
    }
    
    // Eccentricity vector
    double s1 =(1/mu)*(v*v - mu/r);
    double s2 = (1/mu)*r*v_r;
    double ecc_xyz[3] = {s1*Cartesian_vector[0]-s2*Cartesian_vector[3],
                      s1*Cartesian_vector[1]-s2*Cartesian_vector[4],
                      s1*Cartesian_vector[2]-s2*Cartesian_vector[5]};
    // Eccentricity scalar
    double ecc = sqrt(ecc_xyz[0]*ecc_xyz[0] + ecc_xyz[1]*ecc_xyz[1] + ecc_xyz[2]*ecc_xyz[2]);
    // Semi-major axis
    double a = h*h/(mu*(1-ecc*ecc));
    
    // Argument of the pericenter
    double aux =  (-h_xyz[1]*ecc_xyz[0]+h_xyz[0]*ecc_xyz[1])/(N*ecc);
    if (std::fabs(aux) >1.){
        omega = 0.;
    }
    else{
        if (ecc_xyz[2] >=0.){
            omega = std::acos(aux);
        }
        else{
            omega = 2*M_PI - std::acos(aux);
        }
    }
    Orbital_Elements[0]=a;
    Orbital_Elements[1]=ecc;
    Orbital_Elements[2]=i;
    Orbital_Elements[3]=omega;
    Orbital_Elements[4]=Omega;    
}