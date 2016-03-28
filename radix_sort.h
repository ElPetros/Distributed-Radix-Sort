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

#define TEST_CODE 0
#define SEG_FAULT_TEST 0
#define NO_HELPER_FUNC 1    

template <typename T>
std::vector<unsigned int> counting_sort(T* src_begin, T* src_end, T* dst_begin, 
    unsigned int (*key_func)(const T&), unsigned int k = 16, unsigned int d=0);
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
    
    #if SEG_FAULT_TEST
        std::cout << "testing Get Digit " << std::endl;
        std::cout << GET_DIGIT(key_func(*begin), k, 0) << std::endl;
        std::cout << "successful getting digit" << std::endl;
    #endif

    // get comm rank and size
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // The number of elements per processor: n/p
    size_t np = end - begin;
    #if SEG_FAULT_TEST
        std::cout << "in radix sort" << std::endl;
    #endif
    #if TEST_CODE
        std::cout << "SIZE: " << np << "\n";
    #endif

    // the number of histogram buckets = 2^k
    unsigned int num_buckets = 1 << k;

    #if SEG_FAULT_TEST
    std::cout << "in radix sort - before counting sort" << std::endl;
    #endif
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

        int result_index = 0;
        #if SEG_FAULT_TEST
        std::cout << "in radix sort - before histogram" << std::endl;
         #endif
        /* Compute Histogram */
        std::vector<T> result(np);
        std::vector<unsigned int> hist(num_buckets, 0);
        std::vector<unsigned int> L_backup(num_buckets,0);
        
        
    #if NO_HELPER_FUNC
        for (T* iter = begin; iter < end ; iter++){
            hist[GET_DIGIT(key_func(*iter), k, d)]++;
        }

        #if SEG_FAULT_TEST
        std::cout << "in radix sort - after histogram" << std::endl;
         #endif
        std::vector<unsigned int> sum_hist(num_buckets, 0);//Will store the cumulative values
                                                    //Indicates the starting index of each bucket
        for (std::vector<unsigned int>::iterator sum_hist_index = sum_hist.begin() + 1, hist_index = hist.begin() + 1, L_back_index = L_backup.begin() +1; 
            sum_hist_index != sum_hist.end() ; sum_hist_index++, hist_index++, L_back_index++) {
            *sum_hist_index += *(sum_hist_index-1) + *(hist_index-1);
            *L_back_index = *sum_hist_index;

        }
        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after prefix sum" << std::endl;
        #endif
        

        /* Perform Sorting */
        for (T* iter = begin; iter < end ; iter++){//Peforming the sorting
            #if SEG_FAULT_TEST
                std::cout << "sort " << GET_DIGIT(key_func(*iter), k, d) << std::endl;
            #endif
            result[sum_hist[GET_DIGIT(key_func(*iter), k, d)]] = *iter;

            sum_hist[GET_DIGIT(key_func(*iter), k, d)]++;
        }
        #if SEG_FAULT_TEST
                std::cout << "in radix sort - after key based sorting" << std::endl;
         #endif
        
        for (T* iter = begin; iter < end ; iter++, result_index++){//Copying to the original array
            *iter = result[result_index];
        }
        #if SEG_FAULT_TEST
        std::cout << "in radix sort - after copying sorted array" << std::endl;
         #endif
    #else 
        L_backup = counting_sort(begin,end, begin,key_func,k,d);
    #endif

        #if SEG_FAULT_TEST
        std::cout << "in radix sort - before computing parameters to alltoallv" << std::endl;
         #endif

        MPI_Barrier (comm);
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
        
    
        /* Compute G with All Reduce */
        std::vector<unsigned int> G(num_buckets, 0);
        MPI_Allreduce(&hist.front(), &G.front(), hist.size(), MPI_UNSIGNED, MPI_SUM, comm);
        for (std::vector<unsigned int>::iterator iter = G.begin() + 1; iter != G.end() ; iter++){
            *iter += *(iter-1);
        }
        G.insert(G.begin(),0);
        G.pop_back();

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after G" << std::endl;
        #endif

        MPI_Barrier (comm);

        #if TEST_CODE
            std::cout << "G Histogram" << rank << "\n";
            for (std::vector<unsigned int>::iterator iter = G.begin() ; iter != G.end() ; iter++){
                std::cout << *iter << ", ";
            }
            std::cout << "\n";
        #endif
        

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before P" << std::endl;
        #endif
            
        /* Compute P */
        std::vector<unsigned int> P(num_buckets, 0);
        MPI_Exscan(&hist.front(), &P.front(), hist.size(), MPI_UNSIGNED, MPI_SUM, comm);
        MPI_Barrier (comm);
        
        #if TEST_CODE
            std::cout << "P Histogram" << rank << "\n";
            for (std::vector<unsigned int>::iterator iter = P.begin() ; iter != P.end() ; iter++){
                std::cout << *iter << ", ";
            }
            std::cout << "\n";
        #endif

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after P" << std::endl;
        #endif

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before L" << std::endl;
        #endif
        /* Compute L */
        std::vector<unsigned int> L(np, 0);
        int L_index = 0;
        for (T* iter = begin; iter < end ; iter++, L_index++){
            L[L_index] = L_index - L_backup[GET_DIGIT(key_func(*iter), k, d)];
        }
        #if TEST_CODE
            std::cout << "L Histogram" << rank << "\n";
            for (std::vector<unsigned int>::iterator iter = L.begin() ; iter != L.end() ; iter++){
                std::cout << *iter << ", ";
            }
            std::cout << "\n";
        #endif
        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after L" << std::endl;
        #endif

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before T" << std::endl;
        #endif
        /* Compute T */
        std::vector<unsigned int> Tarray(L);
        int T_index = 0;
        for (T* iter = begin; iter < end ; iter++, T_index++){
            Tarray[T_index] += G[GET_DIGIT(key_func(*iter), k, d)] + P[GET_DIGIT(key_func(*iter), k, d)];
        }
        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after T" << std::endl;
        #endif
        MPI_Barrier (comm);
        #if TEST_CODE
                std::cout << "Tarray" << rank << "\n";
                for (std::vector<unsigned int>::iterator iter = Tarray.begin() ; iter != Tarray.end() ; iter++){
                    std::cout << *iter << ", ";
                }
                std::cout << "\n";
        #endif
        

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before Send counts local computation" << std::endl;
        #endif
        /* Compute Send Counts and Send Displacements */
        std::vector<unsigned int> send_counts(p,0);
        std::vector<unsigned int> send_displacements(p,0);
        std::vector<unsigned int> TarrayDividedByLocalSize(Tarray);//For Testing send_counts and send_displacements
        T_index = 1;
        send_counts[*(Tarray.begin())/np]++;
        TarrayDividedByLocalSize[0] /= np;//For testing

        


        for (std::vector<unsigned int>::iterator iter = Tarray.begin()+1 ; iter != Tarray.end() ; iter++, T_index++){
            int destination_processor_rank = (*iter)/np;

            TarrayDividedByLocalSize[T_index] /= np;//For testing

            int previous_destination_processor_rank = *(iter-1)/np;
            send_counts[destination_processor_rank]++;
            if(destination_processor_rank != previous_destination_processor_rank){
                send_displacements[destination_processor_rank] = T_index;
            }
        }

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after Send counts local computation" << std::endl;
        #endif

        #if TEST_CODE
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
        #endif

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before all to all" << std::endl;
        #endif
        //Calculating recv_counts and recv_displacements
        std::vector<unsigned int> recv_counts(p,0);
        MPI_Alltoall(&send_counts.front(), 1, MPI_UNSIGNED, &recv_counts.front(), 1, MPI_UNSIGNED, comm);

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after all to all" << std::endl;
        #endif

        #if TEST_CODE
            // TODO Check if we should be doing all to all with p-1 rather than p.
            // It should be 1. You only receive and send one message per procoess.
            std::cout << "recv_counts" << "\n";
            for (std::vector<unsigned int>::iterator iter = recv_counts.begin() ; iter != recv_counts.end() ; iter++){
                std::cout << *iter << ", ";
            }
            std::cout << "\n";
        #endif

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before recv counts local computation" << std::endl;
        #endif

        /* Compute Receive displacements */
        std::vector<unsigned int> recv_disp(p,0);
        int recv_counter = 0;
        for (std::vector<unsigned int>::iterator iter = recv_disp.begin()+1; iter != recv_disp.end(); iter++) {
            *iter += recv_counts[recv_counter++] + (*(iter-1));
        }

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after recv counts local computation" << std::endl;
        #endif

        #if TEST_CODE
                std::cout << "recv_disp" << "\n";
                for (std::vector<unsigned int>::iterator iter = recv_disp.begin() ; iter != recv_disp.end() ; iter++){
                    std::cout << *iter << ", ";
                }
                std::cout << "\n";
        #endif

        MPI_Barrier(comm);
        // Perform All to AllV
        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before alltoallv" << std::endl;
        #endif
        std::vector<T> recv_buff(begin, end); 
        MPI_Alltoallv(reinterpret_cast<void*>(result.data()), reinterpret_cast<const int*>(send_counts.data()), reinterpret_cast<const int*>(send_displacements.data()), dt, reinterpret_cast<void*>(recv_buff.data()),
                      reinterpret_cast<const int*>(recv_counts.data()), reinterpret_cast<const int*>(recv_disp.data()), dt, comm);

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after alltoallv" << std::endl;
        #endif
        #if TEST_CODE
                std::cout << "recv_buffer" << "\n";


                size_t offset = np*sizeof(T);
                // slightly different line
                for( T* recv_iter = recv_buff.data(); recv_iter < recv_buff.data() + np;  recv_iter++) {
                    std::cout << key_func(*recv_iter) << ", ";
                    // difference in how we obtain the value
                }
                std::cout << std::endl;
        #endif

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - before counting_sort helper" << std::endl;
        #endif
    
        // Locally sort and copy receive buffer using counting sort.
        counting_sort(recv_buff.data(), recv_buff.data() + np, begin, key_func,k,d);

        #if SEG_FAULT_TEST
            std::cout << "in radix sort - after counting_sort helper" << std::endl;
        #endif
        #if TEST_CODE
            std::cout << "\nrecv_buffer SORTED" << "\n";


            offset = np*sizeof(T);
            // slightly different line
            for( T* recv_iter = recv_buff.data(); recv_iter < recv_buff.data() + np;  recv_iter++) {
                std::cout << key_func(*recv_iter) << ", ";
                // difference in how we obtain the value
            }
            std::cout << std::endl;

            std::cout << "Final result\n";
            for(T* iter = begin; iter != end; iter++) {
                
                std::cout << key_func(*iter) << ", ";
            }
            std::cout << std::endl;
        #endif
        MPI_Barrier(comm);
    }
}

template <typename T>
std::vector<unsigned int> counting_sort(T* src_begin, T* src_end, T* dst_begin, unsigned int (*key_func)(const T&), unsigned int k, unsigned int d) {
    size_t np = src_end - src_begin;
    unsigned int num_buckets = 1 << k;
    std::vector<T> result(np);//Temporary sorted vector by the key field
    std::vector<unsigned int> hist(num_buckets, 0);
    std::vector<unsigned int> sum_hist(num_buckets, 0);
    std::vector<unsigned int> L_backup(num_buckets,0);
    std::vector<unsigned int>::iterator sum_hist_index, hist_index;

    for (T* iter = src_begin; iter < src_end ; iter++){//Calculating the histogram
        hist[GET_DIGIT(key_func(*iter), k, d)]++;
    }

    //Indicates the starting index of each bucket
    for (std::vector<unsigned int>::iterator sum_hist_index = sum_hist.begin() + 1, hist_index = hist.begin() + 1, L_back_index = L_backup.begin() +1; 
            sum_hist_index != sum_hist.end() ; sum_hist_index++, hist_index++, L_back_index++) {
            *sum_hist_index += *(sum_hist_index-1) + *(hist_index-1);
            *L_back_index = *sum_hist_index;

        }
    /* Perform Sorting */
    for (T* iter = src_begin; iter < src_end ; iter++){
        result[sum_hist[GET_DIGIT(key_func(*iter), k, d)]] = *iter;
        sum_hist[GET_DIGIT(key_func(*iter), k, d)]++;
    }

    std::copy(result.data(), result.data() + np, dst_begin);    


    return L_backup;
}

