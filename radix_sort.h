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

// returns the value of the digit starting at offset `offset` and containing `k` bits
#define GET_DIGIT(key, k, offset) (((key) >> (offset)) & ((1 << (k)) - 1))

#define TEST_CODE 1
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

    // the number of histogram buckets = 2^k
    unsigned int num_buckets = 1 << k;

    for (unsigned int d = 0; d < 8*sizeof(unsigned int); d += k) {
        std::vector<T> result(np);//Temporary sorted vector by the key field
        std::vector<unsigned int> hist(num_buckets, 0);
        #undef TEST_CODE
        #define TEST_CODE 0
        #if TEST_CODE
            std::cout << k << " " << d << GET_DIGIT(key_func(*begin), k, d) << "\n";
        #endif
        for (T* iter = begin; iter <= end ; iter++){//Calculating the histogram
            hist[GET_DIGIT(key_func(*iter), k, d)]++;
        }

        std::vector<unsigned int> sum_hist(hist);//Will store the cumulative values
        for (std::vector<unsigned int>::iterator next = sum_hist.begin() + 1; next != sum_hist.end() ; next++) {
            *next += *(next-1);
        }

        for (T* iter = begin; iter <= end ; iter++){//Peforming the sorting
            result[sum_hist[key_func(*iter)]-1] = *iter;
        }
        int result_index = 0;
        for (T* iter = begin; iter <= end ; iter++, result_index++){//Copying to the original array
            *iter = result[result_index];
        }

        MPI_Barrier (comm);


        // TODO:
        // 1.) create histogram and sort via bucketing (~ counting sort) DONE
        // 2.) get global histograms (P, G) via MPI_Exscan/MPI_Allreduce,...
        // 3.) calculate send_counts
        // 4.) communicate send_counts to get recv_counts
        // 4.) calculate displacements
        // 6.) MPI_Alltoallv
        // 7.) local sorting via bucketing (~ counting sort)
    }
}

