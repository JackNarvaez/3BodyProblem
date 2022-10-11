/*-------------------------------------------------------------------------------------------------
Reads information about the evolution of the 3-Bodies system and calculates the energy over time.

Author: Narvaez J.
-------------------------------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>  // std::ifstream; // std::ofstream
#include <sstream>  // std::istringstream
#include <stdlib.h>     // atof
#include <cmath>
#include "Integrator.h"

void header(std::ofstream& File, const std::string &name_body, const int &n, const int &k, const int &jump, std::string coord){
    /*---------------------------------------------------------------------------------------------
    Writes the header in file that stores results, as follows:

    name_body
    n k jump

    The name of file is: <coord>_<name_body>.txt
    -----------------------------------------------------------------------------------------------
    Arguments:
    File    :   File where the header is written.
    name_body:  Name of the body.
    n       :   Number of time-iterations in one outer period.
    k       :   Number of outer periods.
    jump    :   Jump size to store data in files.
    coord   :   Keyword.
    ---------------------------------------------------------------------------------------------*/
    File.open("./Files/"+coord+"_"+name_body+".txt");
    File << name_body << std::endl;
    File << n << " " << k << " " << jump << std::endl;
}

void read_evol(const std::string &File_address, double q[][3][6], const int bod, const int it){
    /*---------------------------------------------------------------------------------------------
    Reads information about the evolution of three-body system.
    -----------------------------------------------------------------------------------------------
    Arguments:
    File_address:   File address from which the data is read.
    q   :   Array that stores the evolution of a Body.
    bod :   Number of body.
    it  :   Iterations.
    -----------------------------------------------------------------------------------------------
    Fills the value of inputs as:
        q[i][n] = [xn_i, yn_i, zn_i, vxn_i, vyn_i, vzn_i] for i=1, ..., it. and n=1,2,3.
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

void read_data(const std::string &File_address, const int &n, std::string &name, double & mass){
    /*---------------------------------------------------------------------------------------------
    Reads information about the three-body system.
    -----------------------------------------------------------------------------------------------
    Arguments:
    File_address:   File address from which the data is read.
    n           :   Number of rows in File_address
    name        :   Array that stores the body's name.
    mass        :   mass of the body.
    -----------------------------------------------------------------------------------------------
    Fills the value of inputs as:
        Name = Body's name.
        mass = Mass.
    ---------------------------------------------------------------------------------------------*/
    std::ifstream File;
    File.open (File_address, std::ifstream::in);    // Open file
    std::string line;
    int i = 1;  //Line counter
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
        }
        iss >> data;
        mass = atof(data.c_str());  // Fill mass
        break;
        }
        i++;
    }
    }
}

int main(int argc, char **argv)
{
    std::cout.precision(8);
    std::string name[3];  // Array of names
    name[0] = "Sun";  // Name central body.
    std::string Data = "Data.asc";  // Data file

    int bodies[2] = {atoi(argv[1]), atoi(argv[2])}; // 2nd and 3rd body.

    int k = atoi(argv[3]);    // Number of orbital periods of outer body.
    int jump = atoi(argv[4]); // Stepsize of writing.
    int n = atoi(argv[5]);  // Iterations in a orbital period of outer body.
    double mass[3];// Array of masses.
    mass[0]=1.;
    read_data(Data,bodies[0], name[1], mass[1]);  // Fill data to body 2
    read_data(Data,bodies[1], name[2], mass[2]);  // Fill data to body 
    double q[n][3][6];      // Array of system's evolution.

    std::string Data1 = "./Files/CC_"+name[0]+".txt";  // Data Evolution B1
    std::string Data2 = "./Files/CC_"+name[1]+".txt";  // Data Evolution B2
    std::string Data3 = "./Files/CC_"+name[2]+".txt";  // Data Evolution B3
    double it = k*n/(100*jump);

    read_evol(Data1,q, 0, it);  // Fill B1 data
    read_evol(Data2,q, 1, it);  // Fill B2 data
    read_evol(Data3,q, 2, it);  // Fill B3 data

    std::ofstream E;
    header(E, name[2], n, k, jump, "Energy");

    for (int jj=0; jj<it; jj++){
        E << Energy(q[jj], mass) << std::endl;
    }
    
    E.close();
    return 0;
}