//
// Created by kajetan on 14.06.2020.
//
#include <cmath>
#include <sstream>
#include <limits>
#include <iomanip>
#include "Particle.h"

#define SETPREC std::setprecision(std::numeric_limits<double>::max_digits10)

Vector::Vector(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector::Vector() {
    x = 0.;
    y = 0.;
    z = 0.;
}

Vector Vector::operator+(Vector v) {
    return {x + v.x, y+v.y, z+v.z};
}

Vector Vector::operator-(Vector v) {
    return {x-v.x, y-v.y, z-v.z};
}

double Vector::mod() {
    return std::sqrt(std::pow(x, 2.) + std::pow(y, 2.) + std::pow(z, 2.));
}

Vector Vector::negate() {
    return {-x, -y, -z};
}

Vector Vector::operator*(double v) {
    return {x*v, y*v, z*v};
}

Vector Vector::operator/(double v) {
    return {x/v, y/v, z/v};
}

std::string Vector::toString() {
    std::stringstream stream;
    stream << SETPREC << x << " " << SETPREC << y << " "<< SETPREC << z;
    return stream.str();
}


Particle::Particle(Vector pos, Vector vel) {
    this->pos = pos;
    this->vel = vel;
    this->acc = Vector();
}