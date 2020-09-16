//
// Created by kajetan on 14.06.2020.
//

#include <fstream>
#include <cmath>
#include "Utils.h"
#include<bits/stdc++.h>

#define EPS 4.69041575982343e-08
#define MINIMAL(x) ((std::abs(x) > 1e-10) ? (x) : 1e-10)


double V(double rij, double rik, double rkj){
    double rij_kw, rik_kw, rkj_kw;
    rij_kw = std::pow(rij, 2.);
    rik_kw = std::pow(rik, 2.);
    rkj_kw = std::pow(rkj, 2.);
    double r3 = rij*rik*rkj;
    double lewa = 1./std::pow(r3, 3.);
    double licznik = 3.*(rkj_kw + rik_kw - rij_kw)*(rkj_kw+rij_kw-rik_kw)*(rij_kw+rik_kw-rkj_kw);
    double mianownik = 8.*std::pow(r3, 5.);
    double prawa = licznik/mianownik;
    return lewa + prawa;
}

double V(Vector ri, Vector rj, Vector rk){
    double rij, rik, rkj;
    rij = (rj - ri).mod();
    rij = MINIMAL(rij);
    rik = (rk - ri).mod();
    rik = MINIMAL(rik);
    rkj = (rk - rj).mod();
    rkj = MINIMAL(rkj);

    return V(rij, rik, rkj);

}


Vector dV(Vector ri, Vector rj, Vector rk){
    Vector ret;
    double dVz, z_plus, z_minus;
    volatile double dz;
    double deltaz = (ri.z * EPS);
    z_plus = ri.z + deltaz;
    z_minus = ri.z - deltaz;
    dz = z_plus - z_minus;
    Vector riz_plus(ri.x, ri.y, z_plus);
    Vector riz_minus(ri.x, ri.y, z_minus);
    dVz = V(riz_plus, rj, rk) - V(riz_minus, rj, rk);
    ret.z = std::abs(dz) > EPS*EPS ? dVz/dz : 0;
    double dVy, y_plus, y_minus;
    volatile double dy;
    double deltay = (ri.y * EPS);
    y_plus = ri.y + deltay;
    y_minus = ri.y - deltay;
    dy = y_plus - y_minus;
    Vector riy_plus(ri.x, y_plus, ri.z);
    Vector riy_minus(ri.x, y_minus, ri.z);
    dVy = V(riy_plus, rj, rk) - V(riy_minus, rj, rk);
    ret.y = std::abs(dy) > EPS*EPS ? dVy/dy : 0;
    double dVx, x_plus, x_minus;
    volatile double dx;
    double deltax = (ri.x * EPS) ;
    x_plus = ri.x + deltax;
    x_minus = ri.x - deltax;
    dx = x_plus - x_minus;
    Vector rix_plus(x_plus, ri.y, ri.z);
    Vector rix_minus(x_minus, ri.y, ri.z);
    dVx = V(rix_plus, rj, rk) - V(rix_minus, rj, rk);
    ret.x = std::abs(dx) > EPS*EPS ? dVx/dx : 0;
    return ret * 2.;
}



Vector dV(Particle i, Particle j, Particle k){
    return dV(i.pos, j.pos, k.pos);
}



std::vector<Particle> readFile(const char * filename){
    std::ifstream myfile;
    myfile.open(filename);
    std::vector<Particle> ret;
    double val1, val2, val3;

    while(myfile >> val1){
        myfile >> val2;
        myfile >> val3;
        Vector pos(val1, val2, val3);
        myfile >> val1;
        myfile >> val2;
        myfile >> val3;
        Vector vel(val1, val2, val3);
        ret.emplace_back(pos, vel);
    }
    return ret;
}


Triplet::Triplet(int bid0, int n0,int bid1, int n1,int bid2, int n2) {
    pairs = new std::vector<std::pair<int,int>>();
    this->pairs->reserve(3);
    pairs->push_back({bid0, n0});
    pairs->push_back({bid1, n1});
    pairs->push_back({bid2, n2});
    std::sort(pairs->begin(), pairs->end());
}

bool Triplet::equals(Triplet* other) {
    for(int i=0; i<3; i++){
        if(this->pairs->at(i) != other->pairs->at(i)){
            return false;
        }
    }
    return true;
}

bool Triplet::valid() {
    bool differentBuffer01 = pairs->at(0).first != pairs->at(1).first;
    bool differentBuffer12 = pairs->at(1).first != pairs->at(2).first;
    bool differentIdx01 = pairs->at(0).second != pairs->at(1).second;
    bool differentIdx12 = pairs->at(1).second != pairs->at(2).second;
    bool differentIdx02 = pairs->at(0).second != pairs->at(2).second;
    return (differentBuffer01 && differentBuffer12) || // trzy różne bufory
            (differentBuffer01 && differentIdx12) || // dwa ostatnie takie same
            (differentBuffer12 && differentIdx01) || // dwa pierwsze takie same
            (differentIdx01 && differentIdx12 && differentIdx02); //ten sam bufor
}

TripletBank::TripletBank() {
    this->triplets = std::vector<Triplet*>();
}

void TripletBank::add(int bid0, int n0, int bid1, int n1, int bid2, int n2) {
    triplets.push_back(new Triplet(bid0, n0, bid1, n1, bid2, n2));
}

bool TripletBank::contains(int bid0, int n0, int bid1, int n1, int bid2, int n2) {
    Triplet * t = new Triplet(bid0, n0, bid1, n1, bid2, n2);
    if(!t->valid()){
        return true;
    }
    for(Triplet* tr : triplets){
        if(tr->equals(t)){
            return true;
        }
    }
    return false;
}


