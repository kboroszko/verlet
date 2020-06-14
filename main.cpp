#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>

#define EPS 4.69041575982343e-08
#define MINIMAL(x) (((x) > 1e-10) ? (x) : 1e-10)


double V(double rij, double rik, double rkj){
    double rij_kw, rik_kw, rkj_kw;
    rij_kw = pow(rij, 2.);
    rik_kw = pow(rik, 2.);
    rkj_kw = pow(rkj, 2.);
    double r3 = rij*rik*rkj;
    double lewa = 1./pow(r3, 3.);
    double licznik = 3.*(rkj_kw + rik_kw - rij_kw)*(rkj_kw+rij_kw-rik_kw)*(rij_kw+rik_kw-rkj_kw);
    double mianownik = 8.*pow(r3, 5.);
    double prawa = licznik/mianownik;
    return lewa + prawa;
}

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

Vector::Vector(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector Vector::operator+(Vector v) {
    return {x + v.x, y+v.y, z+v.z};
}

Vector Vector::operator-(Vector v) {
    return {x-v.x, y-v.y, z-v.z};
}

double Vector::mod() {
    return sqrt(pow(x, 2.) + pow(y, 2.) + pow(z, 2.));
}

Vector::Vector() {
    x = 0.;
    y = 0.;
    z = 0.;
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
    stream << std::setprecision(10) << x << " " << std::setprecision(10) << y << " "<< std::setprecision(10)  << z;
    return stream.str();
}



double cos(Vector u, Vector v){
    double licznik = u.x*v.x + u.y*v.y + u.z+v.z;
    double mianownik = u.mod()*v.mod();
    return licznik/mianownik;
}

double Vijk(Vector u, Vector v, Vector w){
    double c1, c2, c3;
    c1 = cos(u, v);
    c2 = cos(w,u);
    c3 = cos(v, w);
    return (1. + 3.*c1*c2*c3)/pow(u.mod()*v.mod()*w.mod(), 3.);
}

double V(Vector ri, Vector rj, Vector rk){
    double rij, rik, rkj;
    rij = (rj - ri).mod();
    rij = MINIMAL(rij);
    rik = (rk - ri).mod();
    rik = MINIMAL(rik);
    rkj = (rk - rj).mod();
    rkj = MINIMAL(rkj);


    double a,b;
    a = V(rij, rik, rkj);
    b = Vijk(rj-ri, rk-ri, rk-rj);

    return a;
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
    ret.z = dz > EPS*EPS ? dVz/dz : 0;
    double dVy, y_plus, y_minus;
    volatile double dy;
    double deltay = (ri.y * EPS);
    y_plus = ri.y + deltay;
    y_minus = ri.y - deltay;
    dy = y_plus - y_minus;
    Vector riy_plus(ri.x, y_plus, ri.z);
    Vector riy_minus(ri.x, y_minus, ri.z);
    dVy = V(riy_plus, rj, rk) - V(riy_minus, rj, rk);
    ret.y = dy > EPS*EPS ? dVy/dy : 0;
    double dVx, x_plus, x_minus;
    volatile double dx;
    double deltax = (ri.x * EPS) ;
    x_plus = ri.x + deltax;
    x_minus = ri.x - deltax;
    dx = x_plus - x_minus;
    Vector rix_plus(x_plus, ri.y, ri.z);
    Vector rix_minus(x_minus, ri.y, ri.z);
    dVx = V(rix_plus, rj, rk) - V(rix_minus, rj, rk);
    ret.x = dx > EPS*EPS ? dVx/dx : 0;
    return ret;
}


class Particle{
public:
    Particle(Vector pos, Vector vel);
    Vector pos;
    Vector vel;
    Vector acc;
};

Particle::Particle(Vector pos, Vector vel){
    this->pos = pos;
    this->vel = vel;
    this->acc = Vector();
}

// 0.00007637217880949 0.00006173286582205 0.0000470984885322


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


Vector dV(Particle i, Particle j, Particle k){
    return dV(i.pos, j.pos, k.pos);
}

//int main(int argc, char *argv[]) {
//    std::ofstream outputFile;
//    outputFile.open ("out.txt");
//    std::ofstream outputFile2;
//    outputFile2.open ("out2.txt");
//
//    Particle p1({-1, 0, 0}, {0, 0, 0});
//    Particle p2({1, 0, 0}, {0, 0, 0});
//    for(int x = -100; x <= 100; x++){
//        for(int y = -100; y <= 100; y++){
//            Particle p3({(double)x/10., (double)y/10., 1}, {0,0,0});
//            outputFile << V(p3.pos, p1.pos, p2.pos) << ((y < 100) ? "," : "");
//        }
//        outputFile << "\n";
//    }
//
//    for(int x = -100; x <= 100; x++){
//        for(int y = -100; y <= 100; y++){
//            Particle p3({(double)x/10., (double)y/10., 0}, {0,0,0});
//            outputFile2 << dV(p3, p1, p2).y << ((y < 100) ? "," : "");
//            if(y < 5 && y > -5){
//                std::cout << x << " " << y << " " << dV(p3, p1, p2).toString() << "\n";
//            }
//        }
//        outputFile2 << "\n";
//    }
//
//    return 0;
//}


//0.000305489 0.000246931 0.000188394
int main(int argc, char *argv[]) {
    if(argc != 5){
        std::cout << "wrong args\n";
        return 1;
    }
    double a = 0.00007637217880949, b = 0.00006173286582205, c = 0.0000470984885322;
//    std::cout << a*4 << " " << b * 4 << " " << c*4 << "\n";
    int n = std::stoi(argv[3]);
    volatile double dt = std::stod(argv[4]);
    auto list = readFile(argv[1]);

    for(Particle p : list){
        std::cout << p.pos.toString() << " " << p.vel.toString() << "\n";
    }

    std::vector<Particle> updated;
    updated = list;

//    //update pos
//    for(int i=0; i<list.size(); i++){
//        list[i].pos = list[i].pos + list[i].vel * dt + list[i].acc * (dt*dt*0.5);
//    }

    for(int i=0; i<list.size(); i++){
//        updated[i].pos = list[i].pos + list[i].vel * dt + list[i].acc * (dt*dt*0.5);
        updated[i].acc = Vector();
        for(int j=0; j<list.size(); j++){
            if( i != j){
                for(int k=j; k<list.size(); k++){
                    if(i != k && j!=k){
                        updated[i].acc = updated[i].acc - dV(list[i], list[j], list[k]);
                    }
                }
            }
        }
//        updated[i].acc = updated[i].acc;
//        updated[i].vel = updated[i].vel + updated[i].acc*dt;
    }

    list = updated;

//update pos
    for(int i=0; i<list.size(); i++){
        list[i].pos = list[i].pos + list[i].vel * dt + list[i].acc * (dt*dt*0.5);
    }

    updated = list;

    for(int i=0; i<list.size(); i++){
//        updated[i].pos = list[i].pos + list[i].vel * dt + list[i].acc * (dt*dt*0.5);
        updated[i].acc = Vector();
        for(int j=0; j<list.size(); j++){
            if( i != j){
                for(int k=j; k<list.size(); k++){
                    if(i != k && j!=k){
                        updated[i].acc = updated[i].acc - dV(list[i], list[j], list[k]);
                    }
                }
            }
        }
//        updated[i].acc = updated[i].acc;
        updated[i].vel = list[i].vel + (list[i].acc + updated[i].acc)*(dt);
    }

    std::ofstream outputFile;
    outputFile.open (argv[2]);
    for(Particle p : updated){
        outputFile << p.pos.toString() << " " << p.vel.toString() << "\n";
        std::cout << p.pos.toString() << " " << p.vel.toString() << "\n";
    }
    outputFile.close();

    return 0;
}















