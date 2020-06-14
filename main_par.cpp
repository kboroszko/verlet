//
// Created by kajetan on 14.06.2020.
//




#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>

#define LOG(msg, rank) {stream<<msg<<"\n";}
#define LOG2(msg1, msg2, rank) {stream<<msg1 << msg2 <<"\n";}
#define LOG3(msg1, msg2, msg3, rank) {stream<<msg1 << msg2 << msg3<<"\n";}


std::stringstream stream;


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

void shift_right(std::vector<double> *bi,int myRank, int numProcesses){
    auto reqs = new MPI_Request[2];
    auto * send_buf = new double[bi->size()];
    auto * recv_buf = bi->data();
    std::copy(recv_buf, recv_buf+bi->size(), send_buf); //duplicate data so async send/rcv possible

    sendRecv(myRank, numProcesses, send_buf, recv_buf, bi->size(), true, reqs);
    MPI_Wait(reqs+1, MPI_STATUS_IGNORE);
}

void shift_left(std::vector<double> *bi,int myRank, int numProcesses){
    auto reqs = new MPI_Request[2];
    auto * send_buf = new double[bi->size()];
    auto * recv_buf = bi->data();
    std::copy(recv_buf, recv_buf+bi->size(), send_buf); //duplicate data so async send/rcv possible

    sendRecv(myRank, numProcesses, send_buf, recv_buf, bi->size(), false, reqs);

    MPI_Waitall(2, reqs, MPI_STATUS_IGNORE);
}


void calculate(std::vector<double>** buffs, int myRank, int numProc){
    int size = (*buffs)[0].size();
    int n0 = (*buffs)[0].at(0);
    int n1 = (*buffs)[1].at(0);
    int n2 = (*buffs)[2].at(0);
    int num_buf0 = n0/size;
    int num_buf1 = n1/size;
    int num_buf2 = n2/size;

    std::cout << num_buf0 << " " << num_buf1 << " " << num_buf2 << "\n";
}

void calculateBufs(std::vector<double>* b0, std::vector<double>* b1, std::vector<double>* b2, int myRank){
    int size = b0->size();
    int n0 = b0->at(0);
    int n1 = b1->at(0);
    int n2 = b2->at(0);
    int num_buf0 = n0/size;
    int num_buf1 = n1/size;
    int num_buf2 = n2/size;

    stream << num_buf0 << " " << num_buf1 << " " << num_buf2 << "\n";
}









int main(int argc, char * argv[])
{
    int myProcessNo;
    int numProcesses;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcessNo);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    //generate buff

    std::vector<double> b0;
    std::vector<double> b1;
    std::vector<double> b2;
    b0.reserve(4);
    b1.reserve(4);
    b2.reserve(4);
    int offset = myProcessNo*4;
    for(int i=0; i<4; i++){
        b0.push_back((double) offset+i);
        b1.push_back((double) offset+i);
        b2.push_back((double) offset+i);
    }

    shift_right(&b0, myProcessNo, numProcesses);
    shift_left(&b2, myProcessNo, numProcesses);


    std::vector<double>* buffs[3];
    buffs[0] = &b0;
    buffs[1] = &b1;
    buffs[2] = &b2;
    int buff_i = 0;

    int p = numProcesses;

    for(int s=p-3; s > 0; s -= 3){
//        LOG2("s=", s, myProcessNo);

        for(int i=0; i < s; i++){
//            LOG2("i=", i, myProcessNo);

            if(i > 0 || s!=p-3){
                shift_right(buffs[buff_i], myProcessNo, numProcesses);
            } else {
                calculateBufs(&b1,&b1,&b1, myProcessNo);
                calculateBufs(&b1,&b1,&b2, myProcessNo);
                calculateBufs(&b0,&b0,&b2, myProcessNo);
            }
            if(s == (p-3)){
                calculateBufs(&b0,&b1,&b1, myProcessNo);
            }

            calculateBufs(&b0,&b1,&b2, myProcessNo);
        }
        buff_i = (buff_i + 1) % 3;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if( (p%3) == 0){
//        LOG("EXTRA", 0)
        shift_right(buffs[buff_i], myProcessNo, numProcesses);
        calculateBufs(&b0,&b1,&b2, myProcessNo);
    }

    for(int i=0; i<numProcesses; i++){
        MPI_Barrier(MPI_COMM_WORLD);
        if(myProcessNo == i){
//            std::cout << "PROC-" << i << "LOG\n";
            std::cout << stream.str() ; //<< "\n";
        }
    }

    MPI_Finalize();

    return 0;
}