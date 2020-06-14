#include <iostream>
#include <fstream>
#include <vector>
#include "Particle.h"
#include "Utils.h"


void updatePos(std::vector<Particle> &list, double dt){
    for(int i=0; i<list.size(); i++){
        list[i].pos = list[i].pos + list[i].vel * dt + list[i].acc * (dt*dt*0.5);
    }
}

void updateVelAndAcc(std::vector<Particle> &list, double dt){
    std::vector<Particle> updated;
    updated = list;


    for(int i=0; i<list.size(); i++){
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
        updated[i].vel = list[i].vel + (list[i].acc + updated[i].acc)*(dt*0.5);
    }

    list = updated;
}



int main(int argc, char *argv[]) {
    if(argc != 5){
        std::cout << "wrong args\n";
        return 1;
    }
    int n = std::stoi(argv[3]);
    volatile double dt = std::stod(argv[4]);
    auto list = readFile(argv[1]);

    std::vector<Particle> updated;
    updated = list;


    for(int i=0; i<list.size(); i++){
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
    }

    list = updated;

    for(int counter = 0; counter < n; counter++){
        updatePos(list, dt);

        updateVelAndAcc(list, dt);
    }


    std::ofstream outputFile;
    outputFile.open (argv[2]);
    for(Particle p : list){
        outputFile << p.pos.toString() << " " << p.vel.toString() << "\n";
        std::cout  << p.pos.toString() << " " << p.vel.toString() << "\n";
    }
    outputFile.close();

    return 0;
}
















