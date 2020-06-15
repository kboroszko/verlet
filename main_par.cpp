//
// Created by kajetan on 14.06.2020.
//




#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <utility>
#include <vector>
#include <sstream>
#include "Particle.h"
#include "Utils.h"

#define LOG(msg, rank) {stream<<msg<<"\n";}
#define LOG2(msg1, msg2, rank) {stream<<msg1 << msg2 <<"\n";}
#define LOG3(msg1, msg2, msg3, rank) {stream<<msg1 << msg2 << msg3<<"\n";}

/**
 * Class to hold particles
 */
class Buffer{
public:
    std::vector<Particle> particles;
    double id;
    Buffer(std::vector<Particle> particles ,double id);
};

Buffer::Buffer(std::vector<Particle> particles, double id) {
    this->particles = std::move(particles);
    this->id = id;
}


std::stringstream stream;
int counter=0;

/**
 * send data from buffFrom to neighbour and recieve from other neighbour buffTo
 * @param myRank
 * @param numProcesses
 * @param buffFrom
 * @param buffTo
 * @param buffSize
 * @param clockwise direction of communication
 * @param reqs
 */
void sendRecv(int myRank, int numProcesses, double* buffFrom, double* buffTo, int buffSize, bool clockwise, MPI_Request* reqs){
    int rankFrom, rankTo;
    if(clockwise){
        rankTo = (myRank + 1) % numProcesses;
        rankFrom = myRank == 0 ? numProcesses - 1 : myRank - 1;
    } else {
        rankFrom = (myRank + 1) % numProcesses;
        rankTo = myRank == 0 ? numProcesses - 1 : myRank - 1;
    }

    MPI_Isend(buffFrom, buffSize, MPI_DOUBLE, rankTo, clockwise, MPI_COMM_WORLD, reqs);
    MPI_Irecv(buffTo, buffSize, MPI_DOUBLE, rankFrom, clockwise, MPI_COMM_WORLD, reqs + 1);
}

/**
 * send buffer to the neighbour on the right
 * @param bi buffer that's content is to be replaced
 * @param myRank
 * @param numProcesses
 */
void shift_right(Buffer *bi,int myRank, int numProcesses){
    stream << myRank << ": sending " << bi->id << " to " << (myRank + 1)%numProcesses << "\n";
    auto reqs = new MPI_Request[2];
    auto * send_buf = new double[bi->particles.size()*9 + 1];
    auto * recv_buf = new double[bi->particles.size()*9 + 1];

    for(int i=0; i<bi->particles.size(); i++){
        int offset = i*9;
        send_buf[offset] = bi->particles.at(i).pos.x;
        send_buf[offset+1] = bi->particles.at(i).pos.y;
        send_buf[offset+2] = bi->particles.at(i).pos.z;
        send_buf[offset+3] = bi->particles.at(i).vel.x;
        send_buf[offset+4] = bi->particles.at(i).vel.y;
        send_buf[offset+5] = bi->particles.at(i).vel.z;
        send_buf[offset+6] = bi->particles.at(i).acc.x;
        send_buf[offset+7] = bi->particles.at(i).acc.y;
        send_buf[offset+8] = bi->particles.at(i).acc.z;
    }
    send_buf[bi->particles.size()*9] = bi->id;

    sendRecv(myRank, numProcesses, send_buf, recv_buf, bi->particles.size()*9+1, true, reqs);


    MPI_Wait(reqs+1, MPI_STATUS_IGNORE);

    for(int i=0; i<bi->particles.size(); i++){
        int offset = i*9;
        bi->particles[i].pos.x = recv_buf[offset];
        bi->particles[i].pos.y = recv_buf[offset+1];
        bi->particles[i].pos.z = recv_buf[offset+2];
        bi->particles[i].vel.x = recv_buf[offset+3];
        bi->particles[i].vel.y = recv_buf[offset+4];
        bi->particles[i].vel.z = recv_buf[offset+5];
        bi->particles[i].acc.x = recv_buf[offset+6];
        bi->particles[i].acc.y = recv_buf[offset+7];
        bi->particles[i].acc.z = recv_buf[offset+8];
    }
    bi->id = recv_buf[bi->particles.size()*9];

}

void shift_left(Buffer *bi,int myRank, int numProcesses){
    auto reqs = new MPI_Request[2];
    auto * send_buf = new double[bi->particles.size()*9 + 1];
    auto * recv_buf = new double[bi->particles.size()*9 + 1];

    for(int i=0; i<bi->particles.size(); i++){
        int offset = i*9;
        send_buf[offset] = bi->particles.at(i).pos.x;
        send_buf[offset+1] = bi->particles.at(i).pos.y;
        send_buf[offset+2] = bi->particles.at(i).pos.z;
        send_buf[offset+3] = bi->particles.at(i).vel.x;
        send_buf[offset+4] = bi->particles.at(i).vel.y;
        send_buf[offset+5] = bi->particles.at(i).vel.z;
        send_buf[offset+6] = bi->particles.at(i).acc.x;
        send_buf[offset+7] = bi->particles.at(i).acc.y;
        send_buf[offset+8] = bi->particles.at(i).acc.z;
    }
    send_buf[bi->particles.size()*9] = bi->id;

    sendRecv(myRank, numProcesses, send_buf, recv_buf, bi->particles.size()*9+1, false, reqs);

    MPI_Wait(reqs+1, MPI_STATUS_IGNORE);

    for(int i=0; i<bi->particles.size(); i++){
        int offset = i*9;
        bi->particles[i].pos.x = recv_buf[offset];
        bi->particles[i].pos.y = recv_buf[offset+1];
        bi->particles[i].pos.z = recv_buf[offset+2];
        bi->particles[i].vel.x = recv_buf[offset+3];
        bi->particles[i].vel.y = recv_buf[offset+4];
        bi->particles[i].vel.z = recv_buf[offset+5];
        bi->particles[i].acc.x = recv_buf[offset+6];
        bi->particles[i].acc.y = recv_buf[offset+7];
        bi->particles[i].acc.z = recv_buf[offset+8];
    }
    bi->id = recv_buf[bi->particles.size()*9];

}

void calculateOne(Buffer *b0){
    int size = b0->particles.size();
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            if(i!=j){
                for(int k=j; k<size; k++){
                    if(i!=k && j!= k){
                        stream << "calc(" << b0->id << ")= " << i << "," << j << "," << k << "\n";
                        b0->particles[i].acc = b0->particles[i].acc - dV(b0->particles[i], b0->particles[j], b0->particles[k]);
                    }
                }
            }
        }
    }
}

void calculateTwo(Buffer *b0, Buffer *b1){
    int size0 = b0->particles.size();
    int size1 = b1->particles.size();
    for(int i=0; i < size0; i++){
        for(int j=0; j<size1; j++){
            for(int k=j+1;k<size1; k++){
                stream << "calc(" << b0->id << "," << b1->id <<")= " << i << "," << j << "," << k << "\n";
                b0->particles[i].acc = b0->particles[i].acc - dV(b0->particles[i], b1->particles[j], b1->particles[k]);
            }
        }
    }
}

void calculateThree(Buffer *b0, Buffer *b1, Buffer *b2){
    int size0 = b0->particles.size();
    int size1 = b1->particles.size();
    int size2 = b2->particles.size();

    if(b0->id == 0 && b1->id == 3 && b2->id == 2){
        stream << "calculateThree :\n";
        stream << "b0[" << b0->id << "]{" << b0->particles[0].pos.toString() << " " << b0->particles[0].vel.toString() << " " << b0->particles[0].acc.toString() << "}\n";
        stream << "b1[" << b1->id << "]{" << b1->particles[0].pos.toString() << " " << b1->particles[0].vel.toString() << " " << b1->particles[0].acc.toString() << "}\n";
        stream << "b2[" << b2->id << "]{" << b2->particles[0].pos.toString() << " " << b2->particles[0].vel.toString() << " " << b2->particles[0].acc.toString()<< "}\n";

    }

    for(int i=0; i<size0; i++){
        for(int j=0; j < size1; j++){
            for(int k=0; k<size2; k++){
                stream << "calc(" << b0->id << "," << b1->id<< "," << b2->id <<")= " << i << "," << j << "," << k << "\n";
                b0->particles[i].acc = b0->particles[i].acc - dV(b0->particles[i], b1->particles[j], b2->particles[k]);
            }
        }
    }
}

/**
 * Method to calculate interactions between particles in three buffers.
 * After this method all particles should have calculated interactions between each other.
 * @param b0
 * @param b1
 * @param b2
 * @param myRank
 */
void calculateBufs(Buffer *b0, Buffer *b1, Buffer* b2, int myRank){
    int n0 = b0->id;
    int n1 = b1->id;
    int n2 = b2->id;
    stream << "calculate " << n0 << " "<< n1 << " " << n2 << '\n';
    if(n0 == n1 && n0 == n2){
        calculateOne(b0);
    } else if(n0 == n1) {
        calculateTwo(b0, b2);
        calculateTwo(b2, b0);
        // calculate n0 n2
        // calculate n1 n2
        // calculate n2 n0
    } else if (n0 == n2){
        calculateTwo(b0, b1);
        calculateTwo(b1, b0);
        // calculate n0 n1
        // calculate n2 n1
        // calculate n1 n0
    } else if (n1 == n2){
        calculateTwo(b0, b1);
        calculateTwo(b1, b0);
        // calculate n0 n1
        // calculate n1 n0
        // calculate n2 n0
    } else {
        calculateThree(b0, b1, b2);
        calculateThree(b1, b0, b2);
        calculateThree(b2, b1, b0);
        // calculate n0 n1 n2
        // calculate n1 n0 n2
        // calculate n2 n0 n1
    }
}

void iterate(Buffer ** buffs, int myProcessNo, int numProcesses){
    int buff_i = 0;
    Buffer *b0 = buffs[0];
    Buffer *b1 = buffs[1];
    Buffer *b2 = buffs[2];

    int p = numProcesses;

    for(int s=p-3; s > 0; s -= 3){
        LOG2("s=", s, myProcessNo);

        for(int i=0; i < s; i++){
            LOG2("i=", i, myProcessNo);

            if(i > 0 || s!=p-3){
                shift_right(buffs[buff_i], myProcessNo, numProcesses);
            } else {
                calculateBufs(b1,b1,b1, myProcessNo);
                calculateBufs(b1,b1,b2, myProcessNo);
                calculateBufs(b0,b0,b2, myProcessNo);
            }
            if(s == (p-3)){
                calculateBufs(b0,b1,b1, myProcessNo);
            }

            calculateBufs(b0,b1,b2, myProcessNo);
        }
        buff_i = (buff_i + 1) % 3;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if( (p%3) == 0){
//        LOG("EXTRA", 0)
        buff_i = buff_i > 0 ? buff_i - 1 : 2;
        shift_right(buffs[buff_i], myProcessNo, numProcesses);

        if((myProcessNo)%(p/3) == 0){
            calculateBufs(b0,b1,b2, myProcessNo);
        }

    }
    //TODO aggregate
}







int main(int argc, char * argv[])
{
    int myProcessNo;
    int numProcesses;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessNo);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    auto list = readFile("test.txt");
    if(list.size() == 0){

        std::cout << myProcessNo << ":error\n" ;

        MPI_Finalize();
        return 1;
    }

    stream << myProcessNo << ":list[i]" << list[myProcessNo].pos.toString()
    << " "<< list[myProcessNo].vel.toString()
    << " "<< list[myProcessNo].acc.toString()
    << "\n";


    std::vector<Particle> updated;
    updated = list;
    int n=1;
    double dt = 0.5;


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

//    list = updated;

    stream << myProcessNo << ":updated" << list.size() << "\n";

    stream << myProcessNo << ":updated[i]" << updated[myProcessNo].pos.toString()
    << " "<< updated[myProcessNo].vel.toString()
    << " "<< updated[myProcessNo].acc.toString()
    << "\n";

    //generate buff
    std::vector<Particle> par;
    par.push_back(list[myProcessNo]);
    Buffer b0(par, myProcessNo);
    Buffer b1(par, myProcessNo);
    Buffer b2(par, myProcessNo);


    stream << myProcessNo << ":" << b0.particles.size()<< b1.particles.size() << b2.particles.size() << "\n";


    stream << "b0[" << b0.id << "]{" << b0.particles[0].pos.toString() << " " << b0.particles[0].vel.toString() << " " << b0.particles[0].acc.toString() << "}\n";
    stream << "b1[" << b1.id << "]{" << b1.particles[0].pos.toString() << " " << b1.particles[0].vel.toString() << " " << b1.particles[0].acc.toString() << "}\n";
    stream << "b2[" << b2.id << "]{" << b2.particles[0].pos.toString() << " " << b2.particles[0].vel.toString() << " " << b2.particles[0].acc.toString()<< "}\n";

    stream << "SHIFT\n";
    shift_right(&b0, myProcessNo, numProcesses);
    shift_left(&b2, myProcessNo, numProcesses);


    Buffer* buffs[3];
    buffs[0] = &b0;
    buffs[1] = &b1;
    buffs[2] = &b2;

    stream << "b0[" << b0.id << "]{" << b0.particles[0].pos.toString() << " " << b0.particles[0].vel.toString() << " " << b0.particles[0].acc.toString() << "}\n";
    stream << "b1[" << b1.id << "]{" << b1.particles[0].pos.toString() << " " << b1.particles[0].vel.toString() << " " << b1.particles[0].acc.toString() << "}\n";
    stream << "b2[" << b2.id << "]{" << b2.particles[0].pos.toString() << " " << b2.particles[0].vel.toString() << " " << b2.particles[0].acc.toString()<< "}\n";

    stream << "ITER\n";

    iterate(buffs, myProcessNo, numProcesses);

    stream << "b0[" << b0.id << "]{" << b0.particles[0].pos.toString() << " " << b0.particles[0].vel.toString() << " " << b0.particles[0].acc.toString() << "}\n";
    stream << "b1[" << b1.id << "]{" << b1.particles[0].pos.toString() << " " << b1.particles[0].vel.toString() << " " << b1.particles[0].acc.toString() << "}\n";
    stream << "b2[" << b2.id << "]{" << b2.particles[0].pos.toString() << " " << b2.particles[0].vel.toString() << " " << b2.particles[0].acc.toString()<< "}\n";

    for(int i=0; i<numProcesses; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        if(myProcessNo == i){
            std::cout << "PROC-" << i << " LOG\n\n";
            std::cout << stream.str() << "\n";
        }
    }

    MPI_Finalize();

    return 0;
}