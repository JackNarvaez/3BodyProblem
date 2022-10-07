#include <iostream>
#include <fstream>  // std::ifstream; // std::ofstream
#include <sstream>  // std::istringstream
#include <stdlib.h>     // atof
#include <cmath>
#include "Integrator.h"

const double G = 4*pow(M_PI,2); // Universal gravitational constant.

double Energy(double q[][6], double mass[]){
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
    rel_vec(q[1],q[0],rel_position[0], 3);
    rel_vec(q[2],q[0],rel_position[1], 3);
    rel_vec(q[2],q[1],rel_position[2], 3);
    for (int ii=0; ii < 3; ii++){
        r[ii] = sqrt(pow(rel_position[ii][0],2)+pow(rel_position[ii][1],2)+pow(rel_position[ii][2],2));
    }
    double U[3] = {mass[0]*mass[1]/r[0], mass[0]*mass[2]/r[1], mass[1]*mass[2]/r[2]};
    double E=0;
    for (int ii=0; ii < 3; ii++){
        E+=0.5*speed2[ii]*mass[ii]-2*G*U[ii];
    }
    return E;
}

void header(std::ofstream& File, const std::string &name_body, const int &n, const int &k, const int &jump, std::string coord){
    /*---------------------------------------------------------------------------------------------
    header:
    Write the header of file to store results, as follows:

    name_body
    n k jump

    The name of file is: coord_name_body.txt
    -----------------------------------------------------------------------------------------------
    Arguments:
        File    :   File to write header.
        name_body:  Name of the body.
        n       :   Number of time-iterations in one outer period.
        k       :   Number of outer periods.
        jump    :   Jump size to store data in files.
        coord   :   OP (Orbital parameters) or CC (Cartesian Coordinates)
    ---------------------------------------------------------------------------------------------*/
    File.open("./Files/"+coord+"_"+name_body+".txt");
    File << name_body << std::endl;
    File << n << " " << k << " " << jump << std::endl;
}

void read_data(const std::string &File_address, double q[][3][6], const int bod, const int it){
    /*---------------------------------------------------------------------------------------------
    read_data:
    Read information about the evolution of three-body system.
    -----------------------------------------------------------------------------------------------
    Arguments:
        File_address    :   File address from which the data is read.
        n               :   Number of rows in File_address
        q               :   Array to store the evolution of a Body.
    -----------------------------------------------------------------------------------------------
    Fill the values inputs as:
        q[i][6] = [xi, yi, zi, vxi, vyi, vzi] for i=1, ..., n
    ---------------------------------------------------------------------------------------------*/
    std::ifstream File;
    File.open (File_address, std::ifstream::in);    // Open file
    std::string line;
    int i = 1;  //Line counter
    std::getline(File,line);
    std::getline(File,line);
	while (!File.eof() && i<it){
	std::getline(File,line);
    std::istringstream iss(line);   // Separate line in columns
    std::string data;
    for(int ii=0; ii < 6; ii++){
        iss >> data;
        q[i][bod][ii] =  atof(data.c_str());   // Fill values of position and velocity
        }
        i++;
    }
}

int main(int argc, char **argv)
{
    std::cout.precision(8);
    std::string name1 = "Sun";  // Name central body.
    std::string name[3] = {"Sun", "Jupiter", "(3040)Kozai"};  // Array of names

    int bodies[2] = {atoi(argv[1]), atoi(argv[2])}; // 2nd and 3rd body.

    int k = atoi(argv[3]);    // Number of orbital periods of outer body.
    int jump = atoi(argv[4]); // Stepsize of writing.
    int n = atoi(argv[5]);  // Iterations in a orbital period of outer body.
    double mass[3];// Array of masses.
    mass[0]=1.;
    mass[1]=9.54792e-04;
    mass[2]=1.0e-16;
    double q[n][3][6];      // Array of system's evolution.

    std::string Data1 = "./Files/CC_"+name[0]+".txt";  // Data Evolution B1
    std::string Data2 = "./Files/CC_"+name[1]+".txt";  // Data Evolution B2
    std::string Data3 = "./Files/CC_"+name[2]+".txt";  // Data Evolution B3
    double it = k*n/(100*jump);

    read_data(Data1,q, 0, it);  // Fill B1 data
    read_data(Data2,q, 1, it);  // Fill B2 data
    read_data(Data3,q, 2, it);  // Fill B3 data

    std::ofstream E;
    header(E, name[2], n, k, jump, "Energy");

    for (int jj=0; jj<it; jj++){
        E << Energy(q[jj], mass) << std::endl;
    }
    
    E.close();
    return 0;
}