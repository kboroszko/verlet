# Verlet algorithm
This is a project for hpc course at UW

Sequential algorithm is located in main_seq.cpp
The parallel version is located in main_par.cpp

The parallel algorithm was described it the publication [here](https://www.researchgate.net/profile/Katherine_Yelick/publication/282380541_A_Computation-_and_Communication-Optimal_Parallel_Direct_3-Body_Algorithm/links/58aaf804aca27206d9bceb90/A-Computation-and-Communication-Optimal-Parallel-Direct-3-Body-Algorithm.pdf).
It generates the unique triplets following algorithm 3 from the paper.

The program assumes your good will, which means that it is not very dilligent in 
checking the validity of the input data and arguments. They just should be correct.

The job for running it on okeanos is in file run.batch. The algorithm assumes that there is at least so many particles as there is processes.