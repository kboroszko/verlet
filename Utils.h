//
// Created by kajetan on 14.06.2020.
//

#ifndef VERLET_UTILS_H
#define VERLET_UTILS_H

#include <vector>
#include "Particle.h"

double V(double rij, double rik, double rkj);
double V(Vector ri, Vector rj, Vector rk);
Vector dV(Vector ri, Vector rj, Vector rk);
Vector dV(Particle i, Particle j, Particle k);

std::vector<Particle> readFile(const char * filename);

#endif //VERLET_UTILS_H
