#include <mpi.h>
#include <vector>







int main(int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	int numbers[4] = {4,3,2,1};	
	int sendcounts[4] = {1,1,1,1};
	int recvcounts[4] = {1,1,1,1};
	int senddisp[4] = {0,1,2,3};
	int recvdisp[4] = {0,1,2,3};
	int recvbuff[4] = {0};
	int rank;
	
	
	std::vector<unsigned int> send_count(sendcounts, sendcounts+4);
	std::vector<unsigned int> recv_count(recvcounts, recvcounts+4);
	std::vector<unsigned int> recv
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	MPI_Alltoallv(numbers, sendcounts, senddisp, MPI_UNSIGNED, recvbuff, recvcounts, recvdisp, MPI_UNSIGNED, MPI_COMM_WORLD);
	
	for(int i = 0; i < 4; i++) {
		std::cout << recvbuff[i] << "( " << rank << "), ";
	}
	
	MPI_Finalize();
	return 0;
}