#ifndef _CONVERSION_H
#define _CONVERSION_H

void op_to_coords(const double &mu, const double &a, const double &ecc, const double &i, const double &omega,
                  const double &Omega, const double &t_0, const double &t, double Cartesian_vector[]);

void coords_to_op(const double &mu, const double Cartesian_vector[], double Orbital_Parameters[]);

#endif