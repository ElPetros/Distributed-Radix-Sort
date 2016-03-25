/**
 * @file    radix_sort.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

/*
 * TODO: implement your radix sort solution in this file
 */

#include <mpi.h>
#include <vector>
#include "mystruct.h"
// returns the value of the digit starting at offset `offset` and containing `k` bits
#define GET_DIGIT(key, k, offset) (((key) >> (offset)) & ((1 << (k)) - 1))

#define TEST_CODE 1
template <typename T>
void counting_sort(T* begin, T* end, unsigned int (*key_func)(const T&), unsigned int k = 16, unsigned int d=0);
/**
 * @brief   Parallel distributed radix sort.
 *
 * This function sorts the distributed input range [begin, end)
 * via lowest-significant-byte-first radix sort.
 *
 * This function will sort elements of type `T`, by the key of type `unsigned int`
 * which is returned by the key-access function `key_func`.
 *
 * The MPI datatype for the templated (general) type `T` has to be passed
 * via the `dt` parameter.
 *
 * @param begin         A pointer to the first element in the local range to be sorted.
 * @param end           A pointer to the end of the range to be sorted. This
 *                      pointer points one past the last element, such that the
 *                      total number of elements is given by `end - begin`.
 * @param key_func      A function with signature: `unsigned int (const T&)`.
 *                      This function returns the key of each element, which is
 *                      used for sorting.
 * @param dt            The MPI_Datatype which represents the type `T`. This
 *                      is used whenever elements of type `T` are communicated
 *                      via MPI.
 * @param comm          The communicator on which the sorting happens.
 *                      NOTE: this is not necessarily MPI_COMM_WORLD. Call
 *                            all MPI functions with this communicator and
 *                            NOT with MPI_COMM_WORLD.
 */
template <typename T>
void radix_sort(T* begin, T* end, unsigned int (*key_func)(const T&), MPI_Datatype dt, MPI_Comm comm, unsigned int k = 16) {
    // get comm rank and size
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // The number of elements per processor: n/p
    size_t np = end - begin;
    std::cout << "SIZE: " << np << "\n";

    // the number of histogram buckets = 2^k
    unsigned int num_buckets = 1 << k;

    for (unsigned int d = 0; d < 8*sizeof(unsigned int); d += k) {
#if TEST_CODE
        std::cout << "*******************Test1" << "\n";
        std::cout << "rank = " << rank << ", k = " << k << ", offset = " << d << "\n";
        std::cout << "UNSorted Local Array" << "\n";
        for (T* iter = begin; iter < end ; iter++){
            std::cout << GET_DIGIT(key_func(*iter), k, d) << ", ";
        }
        std::cout << "\n";
        for (T* iter = begin; iter < end ; iter++){
            std::cout << key_func(*iter) << ", ";
        }
        std::cout << "\n";
#endif
        std::vector<T> result(np);//Temporary sorted vector by the key field
        std::vector<unsigned int> hist(num_buckets, 0);
        for (T* iter = begin; iter < end ; iter++){//Calculating the histogram
            hist[GET_DIGIT(key_func(*iter), k, d)]++;
        }

        std::vector<unsigned int> sum_hist(np, 0);//Will store the cumulative values
                                                    //Indicates the starting index of each bucket
        for (std::vector<unsigned int>::iterator sum_hist_index = sum_hist.begin() + 1, hist_index = hist.begin() + 1; sum_hist_index != sum_hist.end() ; sum_hist_index++, hist_index++) {
            *sum_hist_index += *(sum_hist_index-1) + *(hist_index-1);
        }
        std::vector<unsigned int> L_backup(sum_hist);//sum_hist for calculating L

        for (T* iter = begin; iter < end ; iter++){//Peforming the sorting
            result[sum_hist[GET_DIGIT(key_func(*iter), k, d)]] = *iter;
            sum_hist[GET_DIGIT(key_func(*iter), k, d)]++;
        }
        int result_index = 0;
        for (T* iter = begin; iter < end ; iter++, result_index++){//Copying to the original array
            *iter = result[result_index];
        }
#if TEST_CODE
        std::cout << "*******************Test2" << "\n";
        std::cout << "rank = " << rank << ", k = " << k << ", offset = " << d << "\n";
        std::cout << "Sorted Local Array" << "\n";
        for (T* iter = begin; iter < end ; iter++, result_index++){
            std::cout << GET_DIGIT(key_func(*iter), k, d) << ", ";
        }
        std::cout << "\n";
        for (T* iter = begin; iter < end ; iter++){
            std::cout << key_func(*iter) << ", ";
        }
        std::cout << "\n";
        std::cout << "Histogram" << "\n";
        for (std::vector<unsigned int>::iterator iter = hist.begin() ; iter != hist.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";
#endif
        MPI_Barrier (comm);

        std::vector<unsigned int> G(num_buckets, 0);
        MPI_Allreduce(&hist.front(), &G.front(), hist.size(), MPI_UNSIGNED, MPI_SUM, comm);
        for (std::vector<unsigned int>::iterator iter = G.begin() + 1; iter != G.end() ; iter++){
            *iter += *(iter-1);
        }
        G.insert(G.begin(),0);
        G.pop_back();

#if TEST_CODE
        std::cout << "G Histogram" << rank << "\n";
        for (std::vector<unsigned int>::iterator iter = G.begin() ; iter != G.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";
#endif
        MPI_Barrier (comm);

        std::vector<unsigned int> P(num_buckets, 0);
        MPI_Exscan(&hist.front(), &P.front(), hist.size(), MPI_UNSIGNED, MPI_SUM, comm);
        std::cout << "P Histogram" << rank << "\n";
        for (std::vector<unsigned int>::iterator iter = P.begin() ; iter != P.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";

        MPI_Barrier (comm);

        std::vector<unsigned int> L(np, 0);
        int L_index = 0;
        for (T* iter = begin; iter < end ; iter++, L_index++){
            L[L_index] = L_index - L_backup[GET_DIGIT(key_func(*iter), k, d)];
        }
        std::cout << "L Histogram" << rank << "\n";
        for (std::vector<unsigned int>::iterator iter = L.begin() ; iter != L.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";
        std::vector<unsigned int> Tarray(L);
        int T_index = 0;
        for (T* iter = begin; iter < end ; iter++, T_index++){
            Tarray[T_index] += G[GET_DIGIT(key_func(*iter), k, d)] + P[GET_DIGIT(key_func(*iter), k, d)];
        }
        std::cout << "Tarray" << rank << "\n";
        for (std::vector<unsigned int>::iterator iter = Tarray.begin() ; iter != Tarray.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";
        MPI_Barrier (comm);

        //Calculating send_counts and send_displacements
        std::vector<unsigned int> send_counts(p,0);
        std::vector<unsigned int> send_displacements(p,0);
        std::vector<unsigned int> TarrayDividedByLocalSize(Tarray);//For Testing send_counts and send_displacements
        T_index = 1;
        send_counts[*(Tarray.begin())/np]++;
        TarrayDividedByLocalSize[0] /= np;//For testing
        for (std::vector<unsigned int>::iterator iter = Tarray.begin()+1 ; iter != Tarray.end() ; iter++, T_index++){
            int destination_processor_rank = *iter/np;
            TarrayDividedByLocalSize[T_index] /= np;//For testing
            int previous_destination_processor_rank = *(iter-1)/np;
            send_counts[destination_processor_rank]++;
            if(destination_processor_rank != previous_destination_processor_rank){
                send_displacements[destination_processor_rank] = T_index;
            }
        }
        std::cout << "TarrayDividedByLocalSize" << "\n";
        for (std::vector<unsigned int>::iterator iter = TarrayDividedByLocalSize.begin() ; iter != TarrayDividedByLocalSize.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";
        std::cout << "send_counts" << "\n";
        for (std::vector<unsigned int>::iterator iter = send_counts.begin() ; iter != send_counts.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";
        std::cout << "send_displacements" << "\n";
        for (std::vector<unsigned int>::iterator iter = send_displacements.begin() ; iter != send_displacements.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";

        //Calculating recv_counts and recv_displacements
        std::vector<unsigned int> recv_counts(p,0);
        MPI_Alltoall(&send_counts.front(), p-1, MPI_UNSIGNED, &recv_counts.front(), p-1, MPI_UNSIGNED, comm);
        std::cout << "recv_counts" << "\n";
        for (std::vector<unsigned int>::iterator iter = recv_counts.begin() ; iter != recv_counts.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";


        // Compute Receive displacements
        std::vector<unsigned int> recv_disp(p,0);
        for (std::vector<unsigned int>::iterator iter = recv_disp.begin()+1; iter != recv_disp.end(); iter++) {
            *iter += *(iter-1);
        }

        std::cout << "recv_disp" << "\n";
        for (std::vector<unsigned int>::iterator iter = recv_disp.begin() ; iter != recv_disp.end() ; iter++){
            std::cout << *iter << ", ";
        }
        std::cout << "\n";

        std::vector<T> recv_buff(begin, end); // how to initialize
        MPI_Alltoallv(reinterpret_cast<void*>(result.data()), reinterpret_cast<const int*>(send_counts.data()), reinterpret_cast<const int*>(send_displacements.data()), dt, reinterpret_cast<void*>(recv_buff.data()),
                      reinterpret_cast<const int*>(recv_counts.data()), reinterpret_cast<const int*>(recv_disp.data()), dt, comm);

        std::cout << "recv_buffer" << "\n";


        size_t offset = np*sizeof(T);
        for( T* recv_iter = recv_buff.data(); recv_iter < recv_buff.data() + offset;  recv_iter++) {
            std::cout << key_func(*recv_iter) << ", ";
        }
        std::cout << std::endl;


//        size_t offset = np*sizeof(T);
        // Locally sort receive buffer using counting sort.
//        counting_sort(recv_buff.data(), recv_buff.data() + offset,key_func,k,d);
//
//        // Copy these values to the begin and end
//        int index = 0;
//        T* recv_iter = recv_buff.data();
//        for(T* iter = begin; iter != end; iter++, recv_iter++) {
//            *iter = *recv_iter;
//        }
        // TODO:
        // 1.) create histogram and sort via bucketing (~ counting sort) DONE
        // 2.) get global histograms (P, G) via MPI_Exscan/MPI_Allreduce,... DONE (sorta!)
        // 3.) calculate send_counts locally
        // 4.) communicate send_counts to get recv_counts
        // 4.) calculate displacements
        // 6.) MPI_Alltoallv
        // 7.) local sorting via bucketing (~ counting sort)
    }
}

template <typename T>
void counting_sort(T* begin, T* end, unsigned int (*key_func)(const T&), unsigned int k = 16, unsigned int d=0) {
    size_t np = end - begin;
    unsigned int num_buckets = 1 << k;
    std::vector<T> result(np);//Temporary sorted vector by the key field
    std::vector<unsigned int> hist(num_buckets, 0);
    for (T* iter = begin; iter < end ; iter++){//Calculating the histogram
        hist[GET_DIGIT(key_func(*iter), k, d)]++;
    }

    std::vector<unsigned int> sum_hist(np, 0);//Will store the cumulative values
    //Indicates the starting index of each bucket
    for (std::vector<unsigned int>::iterator sum_hist_index = sum_hist.begin() + 1, hist_index = hist.begin() + 1; sum_hist_index != sum_hist.end() ; sum_hist_index++, hist_index++) {
        *sum_hist_index += *(sum_hist_index-1) + *(hist_index-1);
    }
    std::vector<unsigned int> L_backup(sum_hist);//sum_hist for calculating L

    for (T* iter = begin; iter < end ; iter++){//Peforming the sorting
        result[sum_hist[GET_DIGIT(key_func(*iter), k, d)]] = *iter;
        sum_hist[GET_DIGIT(key_func(*iter), k, d)]++;
    }
    int result_index = 0;

    /* Copy the sorted values to the original array */
    for (T* iter = begin; iter < end ; iter++, result_index++){
        *iter = result[result_index];
    }
}
//void arraySum(void* invec, void* inoutvec, int *len, MPI_Datatype *datatype){
//    MyStruct* in = (MyStruct*) invec;
//    MyStruct* inout = (MyStruct*) inoutvec;
//    for (int i = 0; i < *len; ++i) {
//        inout[i].key += in[i].key;
//    }
//}

