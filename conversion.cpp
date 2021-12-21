// Conversion between orbital parameters and cartesian positions 
// and cartesian velocities.

#include <iostream>
#include <cmath>
#include "conversion.h"

double Halley_method(double x_0, const double &M, const double &ecc, const double &tol, const int NMAX){
    /*-------------------------------------------------------------------------------------------------------
    Halley_method:
    Calculate the root, eccentric anomaly, in the Kepler's equation:
    M = E - e*sin(E),
    by Halley's method.
    ---------------------------------------------------------------------------------------------------------
    Arguments:
        x_0 :   Initial approximation.
        M   :   Mean anomaly.
        ecc :   Eccentricity.
        tol :   Tolerance.
        NMAX:   Maximum number of iterations.
    ---------------------------------------------------------------------------------------------------------
    Return:
        E   :   Eccentric anomaly.
    -------------------------------------------------------------------------------------------------------*/
    double E = std::fmod(x_0, 2*M_PI); // Find the floating point remainder
    double change = 2*tol; //Difference between E_n+1 and E_n [initial value = 2*tol]
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
    /*-------------------------------------------------------------------------------------------------------
    op_to_coords:
    Transform from orbital parameters to cartesian vector.
    Verify the system of units for the arguments.
    ---------------------------------------------------------------------------------------------------------
    Arguments:
        mu  :   G*M, where M = central body's mass.
        a   :   Semi-major axis.
        ecc :   Eccentricity [rad].
        i   :   Inclination [rad].
        omega:  Argument of the pericenter [rad].
        Omega:  Longitude of the ascending node [rad].
        t_0 :   Epoch.
        t   :   Time.
    ---------------------------------------------------------------------------------------------------------
    Change the values of Cartesian_vector to:
        Cartesian_vector[] = {x, y, z, vx, vy, vz}
    -------------------------------------------------------------------------------------------------------*/
    double n = sqrt(mu/pow(a,3));   // Mean motion
    double M = n * (t - t_0);   // Mean anomaly
    double E = Halley_method(0., M, ecc, 10e-10, 100);  // Eccentric anomaly
    double r = a * (1 - ecc * std::cos(E)); // Distance to the central body
    double r_XYZ[3] = {a * (std::cos(E) - ecc), a * sqrt(1 - pow(ecc,2))*std::sin(E), 0.};  // Position vector in orbital frame
    double v_XYZ[3] = {-(sqrt(mu*a)/r)*std::sin(E),(sqrt(mu*a*(1-pow(ecc,2)))/r)*std::cos(E),   
                        0.};    // Velocity vector in the orbital frame
    double R[3][3] = {{std::cos(omega)*std::cos(Omega)-std::sin(omega)*std::sin(Omega)*std::cos(i),
                       -std::sin(omega)*std::cos(Omega)-std::cos(omega)*std::sin(Omega)*std::cos(i),
                       std::sin(Omega)*std::sin(i)},
                      {std::cos(omega)*std::sin(Omega)+std::sin(omega)*std::cos(Omega)*std::cos(i),
                       -std::sin(omega)*std::sin(Omega)+std::cos(omega)*std::cos(Omega)*std::cos(i),
                       -std::cos(Omega)*std::sin(i)},
                      {std::sin(omega)*std::sin(i),std::cos(omega)*std::sin(i),std::cos(i)}};   // Transformation matrix
    // Cartesian vector in inertial frame.
    for (int l=0; l<3; l++){
        Cartesian_vector[l] = R[l][0]*r_XYZ[0]+R[l][1]*r_XYZ[1]+R[l][2]*r_XYZ[2];
        Cartesian_vector[l+3] = R[l][0]*v_XYZ[0]+R[l][1]*v_XYZ[1]+R[l][2]*v_XYZ[2];
    }
}

void coords_to_op(const double &mu, const double Cartesian_vector[], double Orbital_Parameters[]){
    /*-------------------------------------------------------------------------------------------------------
    coords_to_op:
    Transform from Cartesian vector at t=0. to the orbital parameters.
    Verify the system of units for the arguments.
    ---------------------------------------------------------------------------------------------------------
    Arguments: 
        position:  Array with the components of the position in cartesian coordinates [x,y,z].
        Velocity:  Array with the components of the velocity in cartesian coordinates [v_x,v_y,v_z].        
    ---------------------------------------------------------------------------------------------------------
    Change the values of Orbital_Parameters to:
        Orbital_Parameters[i] = xi, where:
        x0 = a:      Semi-major axis.
        x1 = ecc:    Eccentricity [rad].
        x2 = i:      Inclination [rad].
        x3 = omega:  Argument of the pericenter [rad].
        x4 = Omega:  Longitude of the ascending node [rad].
    -------------------------------------------------------------------------------------------------------*/
    // Norm of the position and velocity vectors
    double r = sqrt(pow(Cartesian_vector[0],2)+pow(Cartesian_vector[1],2)+pow(Cartesian_vector[2],2));
    double v = sqrt(pow(Cartesian_vector[3],2)+pow(Cartesian_vector[4],2)+pow(Cartesian_vector[5],2));
    // Radial velocity
    double v_r = (Cartesian_vector[0]*Cartesian_vector[3]+Cartesian_vector[1]*Cartesian_vector[4]+Cartesian_vector[2]*Cartesian_vector[5])/r;
    // Angular momentum
    double h_xyz[3] = {Cartesian_vector[1]*Cartesian_vector[5]-Cartesian_vector[2]*Cartesian_vector[4],
                       Cartesian_vector[2]*Cartesian_vector[3]-Cartesian_vector[0]*Cartesian_vector[5],
                       Cartesian_vector[0]*Cartesian_vector[4]-Cartesian_vector[1]*Cartesian_vector[3]};
    // Norm of the angular momentum
    double h = sqrt(pow(h_xyz[0],2)+pow(h_xyz[1],2)+pow(h_xyz[2],2));
    // Inclination of the orbit
    double i = std::acos(h_xyz[2]/h);
    // Line of Nodes
    double N = sqrt(pow(h_xyz[0],2)+pow(h_xyz[1],2));
    // Longitude of ascending node
    double omega,Omega;
    if (h_xyz[0] >= 0){
        Omega = std::acos(-h_xyz[1]/N);
    }
    else{
        Omega = 2*M_PI - std::acos(-h_xyz[1]/N);
    }
    
    // Eccentricity vector
    double s1 =(1/mu)*(pow(v,2) - mu/r);
    double s2 = (1/mu)*r*v_r;
    double ecc_xyz[3] = {s1*Cartesian_vector[0]-s2*Cartesian_vector[3],
                      s1*Cartesian_vector[1]-s2*Cartesian_vector[4],
                      s1*Cartesian_vector[2]-s2*Cartesian_vector[5]};
    // Eccentricity scalar
    double ecc = sqrt(pow(ecc_xyz[0],2)+pow(ecc_xyz[1],2)+pow(ecc_xyz[2],2));
    // Semi-major axis
    double a = pow(h,2)/(mu*(1-pow(ecc,2)));
    
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
    Orbital_Parameters[0]=a;
    Orbital_Parameters[1]=ecc;
    Orbital_Parameters[2]=i;
    Orbital_Parameters[3]=omega;
    Orbital_Parameters[4]=Omega;    
}