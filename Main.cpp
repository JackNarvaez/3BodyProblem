/*-------------------------------------------------------------------------------------------------
Simulation of a system of 3-bodies interacting with each other gravitationally. 
It finds the time evolution of classical orbital elements (COEs) of the two outer bodies.
It takes by default the Sun and Jupiter as the two massive objects, while the third is a small solar
system object.
The evolution step is performed using the PEFRL integrator, a 4th-order symplectic algorithm.

Author: Narvaez J.
-------------------------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>  // std::ifstream; // std::ofstream
#include <sstream>  // std::istringstream
#include <stdlib.h>     // atof
#include <cmath>
#include <chrono>
#include "conversion.h"
#include "Integrator.h"

void header(std::ofstream& File, const std::string &name_body, const double &dt, const int &n, const int &k, const int &jump, std::string coord);

void read_data(const std::string &File_address, const int &n, std::string &name, double Orbital_Parameters[], double & mass, double &period);

int main(int argc, char **argv)
{
    auto start{std::chrono::steady_clock::now()};
    std::cout.precision(8);
    std::string name1{"Sun"};  // Name central body.
    double body1[6]{0.}; // Cartesian vector of Central body.

    const std::string Data{"Data.asc"};  // Data file
    std::string name[3]{name1};    // Array of names

    int bodies[2]{atoi(argv[1]), atoi(argv[2])}; // 2nd and 3rd body.

    const double G{4*M_PI*M_PI};  // Universal gravitational constant.
    const int k{atoi(argv[3])};    // Number of orbital periods of outer body.
    const int jump{atoi(argv[4])}; // Stepsize of writing.
    const int n{atoi(argv[5])};  // Iterations in a orbital period of outer body.
    double mass[3]{1}, period[3]{0}; // Array of masses and orbital periods, respectively.
    double q[n][3][6]{0};      // Array of system's evolution.
    double body2[7]{0}, body3[7]{0}; // Initial values of orbital elements. 
    double Cartesian_vector2[6]{0}, Cartesian_vector3[6]{0};
    double Orbital_Parameters2[5]{0}, Orbital_Parameters3[5]{0};

    read_data(Data,bodies[0], name[1], body2, mass[1], period[1]);  // Fill data to body 2
    read_data(Data,bodies[1], name[2], body3, mass[2], period[2]);  // Fill data to body 3

    const double T{std::max(period[1],period[2])}; // Orbital period of outer body
    const double dt{T/n}; // Time stepsize
    const double mu{G*mass[0]};

    op_to_coords(mu, body2[0], body2[1], body2[2], body2[3], body2[4], body2[5], body2[6], Cartesian_vector2);
    op_to_coords(mu, body3[0], body3[1], body3[2], body3[3], body3[4], body3[5], body3[6], Cartesian_vector3);
    for (int ll = 0; ll < 6; ll++){
        q[0][0][ll] = body1[ll];
        q[0][1][ll] = Cartesian_vector2[ll];
        q[0][2][ll] = Cartesian_vector3[ll];
    }

    std::ofstream CCBody1, CCBody2, CCBody3, OPBody2, OPBody3;

    header(CCBody1, name[0], dt, n, k, jump, "CC");
    header(CCBody2, name[1], dt, n, k, jump, "CC");
    header(CCBody3, name[2], dt, n, k, jump, "CC");
    header(OPBody2, name[1], dt, n, k, jump, "OP");
    header(OPBody3, name[2], dt, n, k, jump, "OP");
    
    for (int jj=0; jj<k; jj++){
        PEFRL(q[0], mass, n, dt, q);
        for (int ii=0; ii<n; ii+=jump){
            for (int ll=0; ll<6; ll++){
                CCBody1 << q[ii][0][ll] << " ";
                CCBody2 << q[ii][1][ll] << " ";
                CCBody3 << q[ii][2][ll] << " ";
            }
            CCBody1 << '\n';
            CCBody2 << '\n';
            CCBody3 << '\n';
        }
        for (int ee=0; ee<n; ee+=jump){
            rel_vec(q[ee][0], q[ee][1], Cartesian_vector2, 6);
            rel_vec(q[ee][0], q[ee][2], Cartesian_vector3, 6);
            coords_to_op(mu, Cartesian_vector2, Orbital_Parameters2);
            coords_to_op(mu, Cartesian_vector3, Orbital_Parameters3);
            for (int bb = 0; bb < 5; bb++){
                OPBody2 << Orbital_Parameters2[bb] << " ";
                OPBody3 << Orbital_Parameters3[bb] << " ";
            }
            OPBody2 << '\n';
            OPBody3 << '\n';
        }
        for (int ii=0; ii<3; ii++){
            for (int jj=0; jj<6; jj++){
                q[0][ii][jj]=q[n-1][ii][jj];
            }
        }
    }
    
    CCBody1.close();
    CCBody2.close();
    CCBody3.close();
    OPBody2.close();
    OPBody3.close();
    auto end{std::chrono::steady_clock::now()};
    std::cout << "Elapsed time in ms: "
        << std::chrono::duration_cast<std::chrono::seconds>(end-start).count()
        << "\n";
    return 0;
}

void header(std::ofstream& File, const std::string &name_body, const double &dt, const int &n, const int &k, const int &jump, std::string coord){
    /*---------------------------------------------------------------------------------------------
    Writes the header in file that stores results, as follows:

    name_body
    dt n k jump

    The name of file is: <coord>_<name_body>.txt
    -----------------------------------------------------------------------------------------------
    Arguments:
    File    :   File where the header is written.
    name_body:  Name of the body.
    dt      :   Discrete time step.
    n       :   Number of time-iterations in one outer period.
    k       :   Number of outer periods.
    jump    :   Jump size to store data in files.
    coord   :   OP (Orbital parameters) or CC (Cartesian Coordinates)
    ---------------------------------------------------------------------------------------------*/
    File.open("./Files/"+coord+"_"+name_body+".txt");
    File << name_body << std::endl;
    File << dt << " " << n << " " << k << " " << jump << std::endl;
}

void read_data(const std::string &File_address, const int &n, std::string &name, double Orbital_Parameters[], double & mass, double &period){
    /*---------------------------------------------------------------------------------------------
    Reads information about the three-body system.
    -----------------------------------------------------------------------------------------------
    Arguments:
    File_address    :   File address from which the data is read.
    n               :   Number of rows in File_address
    name            :   Array that stores the body's name.
    Orbital_Parameters: Array that stores orbital elements.
    mass            :   mass of the body.
    period          :   Orbital period.
    -----------------------------------------------------------------------------------------------
    Fills the value of inputs as:
        Name = Body's name.
        Orbital_Parameters[] = {a, ecc, i, omega, Omega, epoch, t}
        mass = Mass.
        period = Orbital period.
    ---------------------------------------------------------------------------------------------*/
    std::ifstream File;
    File.open (File_address, std::ifstream::in);    // Open file
    std::string line;
    int i{1};  //Line counter
	while (!File.eof()){
	std::getline(File,line);
    // Omit empty lines and comments
	if (line.length() == 0 || line[0] == '#'){
		continue;
    }else{
        if(i==n){ // Find line of body in data file.
        std::istringstream iss(line);   // Separate line in columns
        std::string data;
        iss >> name;    // Fill name
        for(int ii=0; ii < 7; ii++){
            iss >> data;
            Orbital_Parameters[ii] =  atof(data.c_str());   // Fill Orbital elements 
        }
        iss >> data;
        mass = atof(data.c_str());  // Fill mass
        iss >> data;
        period = atof(data.c_str());    // Fill orbital period
        break;
        }
        i++;
    }
    }
}