#include <interface.hpp>
#include <iostream>

int main(int argc, char** argv){
    int nworkers, startloop, endloop, step, chunksize;

    if (argc != 6) {
        std::cerr << "use: "
                  << argv[0]
                  << " nworkers startloop endloop step chunksize\n";
        return -1;
    }
    nworkers = atoi(argv[1]);
    startloop = atoi(argv[2]);
    endloop = atoi(argv[3]);
    step = atoi(argv[4]);
    chunksize = atoi(argv[5]);


    std::vector<uint> v;
    nornir::ParallelFor pf(nworkers, new nornir::Parameters("parameters.xml"));
    for(size_t i = 0; i < 10; i++){
        v.resize(endloop - startloop);
        std::cout << "outer iteration " << i << std::endl;
        pf.parallel_for(startloop, endloop, step, chunksize,
        [&v](long long int idx, long long int id){
           v[idx] = idx;
        });

        for(uint x : v){
            std::cout << x << std::endl;
        }
        v.clear();
    }
    return 0;
}
