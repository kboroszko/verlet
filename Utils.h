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

class Triplet{
public:
    std::vector<std::pair<int,int>>* pairs;
    Triplet(int bid0, int n0,int bid1, int n1,int bid2, int n2);
    bool equals(Triplet* other);
    bool valid();
};

class TripletBank {
public:
    std::vector<Triplet*> triplets;
    TripletBank();
    void add(int bid0, int n0,int bid1, int n1,int bid2, int n2);
    bool contains(int bid0, int n0,int bid1, int n1,int bid2, int n2);
};

#endif //VERLET_UTILS_H
