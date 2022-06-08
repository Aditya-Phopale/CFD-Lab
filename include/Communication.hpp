#include<mpi.h>
class Communication{
    public: 
        static void init_parallel(int argn, char **args, int rank, int size);
        static void finalize();
        static void communicate();
        static double reduce_min();
        static double reduce_sum();
};