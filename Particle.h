//
// Created by kajetan on 14.06.2020.
//

#ifndef VERLET_PARTICLE_H
#define VERLET_PARTICLE_H


#include <string>

class Vector{
public:
    Vector();
    Vector(double x, double y, double z);
    Vector operator +(Vector v);
    Vector operator -(Vector v);
    Vector operator /(double v);
    Vector operator *(double v);
    double mod();
    Vector negate();
    std::string toString();
    volatile double x, y, z;
};


class Particle{
public:
    Particle(Vector pos, Vector vel);
    Particle();
    Vector pos;
    Vector vel;
    Vector acc;
};


#endif //VERLET_PARTICLE_H
