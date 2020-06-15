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
#include <zconf.h>
#include <fstream>
#include "Particle.h"
#include "Utils.h"

#define LOG(msg, rank) {std::cout<<msg<<"\n";}
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


//std::stringstream stream;

int numParticles = 0;
int numProcesses = 0;

int getFirstParticleIndexOfProcess(int myRank) {
    if(numParticles == 0 || numProcesses == 0){
        throw std::runtime_error("numParticles or numProceses not set!");
    }
    if(myRank == numProcesses){
        return numParticles;
    }
    int rowsInEveryone = numParticles / numProcesses;
    int rest = numParticles % rowsInEveryone;

    if(rest == 0){
        return myRank * rowsInEveryone;
    } else {
        if(myRank <= rest){
            return myRank * (rowsInEveryone + 1);
        } else {
            int fullRows = rest * (rowsInEveryone + 1);
            return fullRows + (myRank - rest) * rowsInEveryone;
        }
    }
}

int getBufSizeOfProcess(int myRank){
    if(numParticles == 0 || numProcesses == 0){
        throw std::runtime_error("numParticles or numProceses not set!");
    }
    int a = getFirstParticleIndexOfProcess(myRank+1);
    int b = getFirstParticleIndexOfProcess(myRank);
    return a - b;
}

int getMaxBufSize(){
    if(numParticles == 0 || numProcesses == 0){
        throw std::runtime_error("numParticles or numProceses not set!");
    }
    return (numParticles + numProcesses - 1) / numProcesses;
}

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
void shift_right(Buffer *bi,int myRank){

    if(numParticles == 0 || numProcesses == 0){
        throw std::runtime_error("numParticles or numProceses not set!");
    }


//    stream << myRank << ": sending " << bi->id << " to " << (myRank + 1)%numProcesses << "\n";
    auto reqs = new MPI_Request[2];
    int maxSize = getMaxBufSize();
    auto * send_buf = new double[maxSize*9 + 1];
    auto * recv_buf = new double[maxSize*9 + 1];


    for(int i=0; i<getBufSizeOfProcess(bi->id); i++){
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
    send_buf[maxSize*9] = bi->id;

    sendRecv(myRank, numProcesses, send_buf, recv_buf, maxSize*9+1, true, reqs);


    MPI_Wait(reqs+1, MPI_STATUS_IGNORE);

    bi->id = recv_buf[maxSize*9];

    for(int i=0; i<getBufSizeOfProcess(bi->id); i++){
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

}

void shift_left(Buffer *bi,int myRank){

    if(numParticles == 0 || numProcesses == 0){
        throw std::runtime_error("numParticles or numProceses not set!");
    }

    int maxSize = getMaxBufSize();
    auto reqs = new MPI_Request[2];
    auto * send_buf = new double[maxSize*9 + 1];
    auto * recv_buf = new double[maxSize*9 + 1];

    for(int i=0; i<getBufSizeOfProcess(bi->id); i++){
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
    send_buf[maxSize*9] = bi->id;

    sendRecv(myRank, numProcesses, send_buf, recv_buf, maxSize*9+1, false, reqs);

    MPI_Wait(reqs+1, MPI_STATUS_IGNORE);


    bi->id = recv_buf[maxSize*9];

    for(int i=0; i<getBufSizeOfProcess(bi->id); i++){
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

}

/**
 * calculate interactions between particles in buff b0
 * @param b0
 */
void calculateOne(Buffer *b0){
    int size = getBufSizeOfProcess(b0->id);
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            if(i!=j){
                for(int k=j; k<size; k++){
                    if(i!=k && j!= k){
//                        stream << "calc(" << b0->id << ")= " << i << "," << j << "," << k << "\n";
                        b0->particles[i].acc = b0->particles[i].acc - dV(b0->particles[i], b0->particles[j], b0->particles[k]);
                    }
                }
            }
        }
    }
}

/**
 * Calculate interactions between particles in buffers b0 and b1.
 * Note that only particles in buffer b0 will have updated acc vector.
 * @param b0
 * @param b1
 */
void calculateTwo(Buffer *b0, Buffer *b1){
    int size0 = getBufSizeOfProcess(b0->id);
    int size1 = getBufSizeOfProcess(b1->id);
    for(int i=0; i < size0; i++){
        for(int j=0; j<size1; j++){
            for(int k=j+1;k<size1; k++){
//                stream << "calc(" << b0->id << "," << b1->id <<")= " << i << "," << j << "," << k << "\n";
                b0->particles[i].acc = b0->particles[i].acc - dV(b0->particles[i], b1->particles[j], b1->particles[k]);
            }
        }
    }
}

/**
 * Calculate interactions between particles in buffers b0 and b1.
 * Note that only particles in buffer b0 will have updated acc vector.
 * @param b0
 * @param b1
 * @param b2
 */
void calculateThree(Buffer *b0, Buffer *b1, Buffer *b2){
    int size0 = getBufSizeOfProcess(b0->id);
    int size1 = getBufSizeOfProcess(b1->id);
    int size2 = getBufSizeOfProcess(b2->id);

    for(int i=0; i<size0; i++){
        for(int j=0; j < size1; j++){
            for(int k=0; k<size2; k++){
//                stream << "calc(" << b0->id << "," << b1->id<< "," << b2->id <<")= " << i << "," << j << "," << k << "\n";
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
void calculateBufs(Buffer *b0, Buffer *b1, Buffer* b2){
    int n0 = b0->id;
    int n1 = b1->id;
    int n2 = b2->id;
//    stream << "calculate " << n0 << " "<< n1 << " " << n2 << '\n';
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
//        LOG2("s=", s, myProcessNo);

        for(int i=0; i < s; i++){
//            LOG2("i=", i, myProcessNo);

            if(i > 0 || s!=p-3){
                shift_right(buffs[buff_i], myProcessNo);
            } else {
                calculateBufs(b1,b1,b1);
                calculateBufs(b1,b1,b2);
                calculateBufs(b0,b0,b2);
            }
            if(s == (p-3)){
                calculateBufs(b0,b1,b1);
            }

            calculateBufs(b0,b1,b2);
        }
        buff_i = (buff_i + 1) % 3;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if( (p%3) == 0){
//        LOG("EXTRA", 0)
        buff_i = buff_i > 0 ? buff_i - 1 : 2;
        shift_right(buffs[buff_i], myProcessNo);

        if((myProcessNo)%(p/3) == 0){
            calculateBufs(b0,b1,b2);
        }

    }

}

double * aggregate(Buffer ** buffs, int myProcessNo){

    Buffer *b0 = buffs[0];
    Buffer *b1 = buffs[1];
    Buffer *b2 = buffs[2];

    if( b0->id > b1->id ){
        Buffer * tmp;
        tmp = b1;
        b1 = b0;
        b0 = tmp;
    }
    if( b2->id < b1->id ){
        Buffer * tmp;
        tmp = b1;
        b1 = b2;
        b2 = tmp;
    }
    if( b0->id > b1->id ){
        Buffer * tmp;
        tmp = b1;
        b1 = b0;
        b0 = tmp;
    }

    Buffer * ordered[3];
    ordered[0] = b0;
    ordered[1] = b1;
    ordered[2] = b2;

    auto * aggr_acc = new double [3*getBufSizeOfProcess(myProcessNo)];

    for(int c=0; c<numProcesses; c++){
        double send_acc[3 * getMaxBufSize()];

        if(c == b0->id || c == b1->id || c == b2->id){
            Buffer* bi;
            if(c == b0->id ){
                bi = b0;
            } else if(c == b1->id){
                bi = b1;
            } else {
                bi = b2;
            }
//            bi->particles.resize(getMaxBufSize());
            for(int i=0; i<getBufSizeOfProcess(c); i++){
                int offset = i*3;
                send_acc[offset] = bi->particles.at(i).acc.x;
                send_acc[offset + 1] = bi->particles.at(i).acc.y;
                send_acc[offset + 2] = bi->particles.at(i).acc.z;
            }
        } else {
            for(int i=0; i < 3 * getBufSizeOfProcess(c); i++){
                send_acc[i] = 0.;
            }
        }

        if(myProcessNo == c){
            MPI_Reduce(send_acc, aggr_acc, 3 * getBufSizeOfProcess(c),
                       MPI_DOUBLE, MPI_SUM, myProcessNo,
                       MPI_COMM_WORLD);
        } else {
            MPI_Reduce(send_acc, nullptr, 3 * getBufSizeOfProcess(c),
                       MPI_DOUBLE, MPI_SUM, c,
                       MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    return aggr_acc;
}

void update(Buffer * myBuffer, int myProcessNo, double dt){
    std::vector<Particle> oldList;
    oldList = myBuffer->particles;
    int size = getBufSizeOfProcess(myProcessNo);
    for(int i=0; i < size; i++){
        myBuffer->particles[i].pos = myBuffer->particles[i].pos + myBuffer->particles[i].vel * dt + myBuffer->particles[i].acc * (dt*dt*0.5);
        myBuffer->particles[i].acc = Vector();
    }

    Buffer b0(myBuffer->particles, myProcessNo);
    b0.particles.resize(getMaxBufSize());
    Buffer b1(myBuffer->particles, myProcessNo);
    b1.particles.resize(getMaxBufSize());
    Buffer b2(myBuffer->particles, myProcessNo);
    b2.particles.resize(getMaxBufSize());

    shift_right(&b0, myProcessNo);
    shift_left(&b2, myProcessNo);

    Buffer* buffs[3];
    buffs[0] = &b0;
    buffs[1] = &b1;
    buffs[2] = &b2;

    iterate(buffs, myProcessNo, numProcesses);

    MPI_Barrier(MPI_COMM_WORLD);
    double * aggr_acc = aggregate(buffs, myProcessNo);

    Buffer * bi = myBuffer;

    for(int i=0; i<size; i++){
        int offset = i*3;
        bi->id = myProcessNo;
        bi->particles.at(i).vel.x = oldList.at(i).vel.x +(aggr_acc[offset] + oldList.at(i).acc.x)*(dt*0.5);
        bi->particles.at(i).vel.y = oldList.at(i).vel.y +(aggr_acc[offset + 1] + oldList.at(i).acc.y)*(dt*0.5);
        bi->particles.at(i).vel.z = oldList.at(i).vel.z +(aggr_acc[offset + 2]+ oldList.at(i).acc.z)*(dt*0.5);
        bi->particles.at(i).acc.x = aggr_acc[offset];
        bi->particles.at(i).acc.y = aggr_acc[offset + 1];
        bi->particles.at(i).acc.z = aggr_acc[offset + 2];
    }
}

void updateAccFirst(Buffer * myBuffer, int myProcessNo){
    Buffer b0(myBuffer->particles, myProcessNo);
    b0.particles.resize(getMaxBufSize());
    Buffer b1(myBuffer->particles, myProcessNo);
    b1.particles.resize(getMaxBufSize());
    Buffer b2(myBuffer->particles, myProcessNo);
    b2.particles.resize(getMaxBufSize());

    shift_right(&b0, myProcessNo);
    shift_left(&b2, myProcessNo);

    Buffer* buffs[3];
    buffs[0] = &b0;
    buffs[1] = &b1;
    buffs[2] = &b2;

    iterate(buffs, myProcessNo, numProcesses);

    MPI_Barrier(MPI_COMM_WORLD);
    double * aggr_acc = aggregate(buffs, myProcessNo);

    Buffer * bi = myBuffer;

    for(int i=0; i<getBufSizeOfProcess(myProcessNo); i++){
        int offset = i*3;
        bi->id = myProcessNo;
        bi->particles.at(i).acc.x = aggr_acc[offset];
        bi->particles.at(i).acc.y = aggr_acc[offset + 1];
        bi->particles.at(i).acc.z = aggr_acc[offset + 2];
    }
}

void printToFile(Buffer * buff, int myRank, int iter, std::string filename){
    std::ofstream outputFile;
    filename += "_";
    filename += std::to_string(iter);
    filename += ".txt";
    for(int i=0; i < numProcesses; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        if(myRank == i){

            if(myRank == 0){
                outputFile.open(filename);
            } else {
                outputFile.open(filename, std::ios_base::app);
            }

            for(Particle p : buff->particles){
                outputFile << p.pos.toString() << " " << p.vel.toString() << "\n";
            }

            outputFile.close();
        }
    }
}

//void printBuff(Buffer *b){
//    stream << "{" ;
//    for(int i=0; i<getBufSizeOfProcess(b->id); i++){
//        stream << "{" << b->particles[i].pos.toString() << " "
//        << b->particles[i].vel.toString() << " "
//        << b->particles[i].acc.toString() << "}";
//    }
//    stream << "}";
//}

int main(int argc, char * argv[])
{
    if(argc != 5 && argc != 6){
        std::cout << "wrong args \n";
        return 1;
    }

    bool verbose = argc == 6;

    int myProcessNo;

    int n = std::stoi(argv[3]);
    volatile double dt = std::stod(argv[4]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessNo);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

//    if(myProcessNo == 0){
//        int xx = 0;
//        while(xx == 0)
//            sleep(1);
//    }

    std::vector<Particle> list;

    if(myProcessNo == 0){
        list = readFile(argv[1]);
        if(list.empty()){
            throw std::runtime_error("FILE EMPTY");
        }

        numParticles = list.size();
    }


    MPI_Bcast(
            &numParticles,
            1,
            MPI_INT,
            0,
            MPI_COMM_WORLD);

    int maxSize = getMaxBufSize();


    std::vector<Particle> par;

    if(myProcessNo == 0){

        for(int i=0; i<getBufSizeOfProcess(myProcessNo); i++){
            par.push_back(list[i]);
        }

        for(int c=1; c<numProcesses; c++){
            auto * send_buf = new double[maxSize*6];

            for(int i=getFirstParticleIndexOfProcess(c); i<getFirstParticleIndexOfProcess(c+1); i++){
                int offset = (i-getFirstParticleIndexOfProcess(c))*6;
                send_buf[offset] = list.at(i).pos.x;
                send_buf[offset+1] = list.at(i).pos.y;
                send_buf[offset+2] = list.at(i).pos.z;
                send_buf[offset+3] = list.at(i).vel.x;
                send_buf[offset+4] = list.at(i).vel.y;
                send_buf[offset+5] = list.at(i).vel.z;
            }

            MPI_Send(
                    send_buf,
                    maxSize*6,
                    MPI_DOUBLE,
                    c,
                    0,
                    MPI_COMM_WORLD);
        }
    } else {
        //recieve data from proc 0
        auto * recv_buf = new double[maxSize*6];

        MPI_Recv(
                recv_buf,
                maxSize*6,
                MPI_DOUBLE,
                0,
                0,
                MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);


        for(int i=0; i<getBufSizeOfProcess(myProcessNo); i++){
            int offset = i*6;
            Vector pos = Vector(recv_buf[offset], recv_buf[offset+1], recv_buf[offset+2]);
            Vector vel = Vector(recv_buf[offset+3], recv_buf[offset+4], recv_buf[offset+5]);
            par.emplace_back(pos, vel);
        }
    }

    par.resize(getMaxBufSize());
    auto * myBuff = new Buffer(par, myProcessNo);
//    stream << "Before: myBuff[" << myBuff->id << "]";
//    printBuff(myBuff);
//    stream << "\n";

//    LOG("OK", myProcessNo)

    updateAccFirst(myBuff, myProcessNo);
    for(int i=0; i < n; i++){
        update(myBuff, myProcessNo, dt);

        if(verbose && i < n-1){
            printToFile(myBuff, myProcessNo, i+1,argv[2]);
        }
    }
    printToFile(myBuff, myProcessNo, n, argv[2]);


//    stream << "After: myBuff[" << myBuff->id << "]";
//    printBuff(myBuff);
//    stream << "\n";

//    for(int i=0; i<numProcesses; i++){
//        MPI_Barrier(MPI_COMM_WORLD);
//        if(myProcessNo == i){
//            std::cout << "PROC-" << i << " LOG\n\n";
//            std::cout << stream.str() << "\n";
//        }
//    }


    MPI_Finalize();

    return 0;
}