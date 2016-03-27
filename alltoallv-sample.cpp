#include <mpi.h>
#include <vector>







int main(int argc, char* argv[]) {
	MPI_Init(NULL, NULL);

	int size = 5;
	int numbers[5] = {4,3,2,1,0};	
	int sendcounts[5] = {1,1,1,1,1};
	int recvcounts[5] = {1,1,1,1,1};
	int senddisp[5] = {0,1,2,3,4};
	int recvdisp[5] = {0,1,2,3,4};
	int recvbuff[5] = {0};
	int rank;
	
	
	std::vector<unsigned int> send_count(sendcounts, sendcounts+size);
	std::vector<unsigned int> recv_count(recvcounts, recvcounts+size);
	std::vector<unsigned int> recv_disp(recvdisp, recvdisp+size);
	std::vector<unsigned int> send_disp(senddisp, senddisp+size);
	std::vector<unsigned int> numbers_vec(numbers, numbers+size);
	std::vector<unsigned int> recv_vec(numbers_vec);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Alltoallv(reinterpret_cast<void*>(numbers_vec.data()), reinterpret_cast<const int*>(send_count.data()), 
		reinterpret_cast<const int*>(send_disp.data()), MPI_UNSIGNED, reinterpret_cast<void*>(recv_vec.data()), reinterpret_cast<const int*>(recv_count.data()), reinterpret_cast<const int*>(recv_disp.data()), MPI_UNSIGNED, MPI_COMM_WORLD);
	
	
	size_t offset = 1*sizeof(unsigned int);
	for(unsigned int* recv_iter = recv_vec.data(); recv_iter < recv_vec.data( ) + offset; recv_iter++) {
		std::cout << *recv_iter << "( " << rank << "), ";
	}
	std::cout << std::endl;
	
	MPI_Finalize();
	return 0;
}