#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define	NZ	3
#define NX	4
#define NT	2
#define MASTER	0
#define TAG 100

int main(int argc, char **argv)
{
	int numProcesses, taskId;
	int i, j, k, irank;
	char global[NT*NZ*NX] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X'};
	char *local;
	int avg_face, blk_face, rem_face;
	int *sendcounts;
	int *displs;
	int offset;
	int recvcount;
	MPI_Datatype type_ztPlane, type_ztPlane_resized, type_ztxBlock, type_ztxBlock_resized;
	MPI_Status status;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	if (taskId == 0)
	{
		printf("Printing global\n");
		for (k = 0; k < NT; k++)
		{
			for (i = 0; i < NZ; i++)//j = 0; j < NX; j++
			{
				for (j = 0; j < NX; j++)//i = 0; i < NZ; i++
				{
					printf("global[%d] = %c\t", k*(NZ*NX)+j*NZ+i, global[k*(NZ*NX)+j*NZ+i]);
				}
				printf("\n");
			}
			printf("\n\n");
		}
	}
	avg_face = NX / numProcesses;
	rem_face = NX % numProcesses;
	offset = 0;
	sendcounts = (int*)calloc(numProcesses, sizeof(int));
	displs = (int*)calloc(numProcesses, sizeof(int));
	for (irank = 0; irank < numProcesses; irank++)
	{
		blk_face = (irank < rem_face) ? (avg_face + 1) : (avg_face);
		sendcounts[irank] = blk_face;
		displs[irank] = offset;
		offset += sendcounts[irank];
	}
	recvcount = sendcounts[taskId];

	
	if (taskId == 0)
	{
		// send datatype
		MPI_Type_vector(NT, NZ, NZ*NX, MPI_CHAR, &type_ztPlane);
		MPI_Type_commit(&type_ztPlane);
		MPI_Type_create_resized(type_ztPlane, 0, NZ * sizeof(char), &type_ztPlane_resized);
		MPI_Type_commit(&type_ztPlane_resized);
	}
    
    // receive datatype
    MPI_Type_vector(NT, NZ, NZ*recvcount, MPI_CHAR, &type_ztxBlock);
    MPI_Type_commit(&type_ztxBlock);
    MPI_Type_create_resized(type_ztxBlock, 0, NZ * sizeof(char), &type_ztxBlock_resized);
    MPI_Type_commit(&type_ztxBlock_resized);

	printf("rank:%d, sendcounts[%d] = %d, displs[%d] = %d, recvcount = %d\n", 
		taskId, taskId, sendcounts[taskId], taskId, displs[taskId], recvcount);
	local = (char*)calloc(recvcount * NZ * NT, sizeof(char));
    // MPI_Scatterv(global, sendcounts, displs, type_ztPlane_resized,
    //         local, recvcount * NZ * NT, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    MPI_Scatterv(global, sendcounts, displs, type_ztPlane_resized,
            local, recvcount, type_ztxBlock_resized, MASTER, MPI_COMM_WORLD);
 	printf("Printing local for rank %d\n", taskId);
	for (k = 0; k < NT; k++)
	{
		for (i = 0; i < NZ; i++)//j = 0; j < NX; j++
		{
			for (j = 0; j < recvcount; j++)//i = 0; i < NZ; i++
			{
				printf("local[%d] = %c\t", k*(NZ*recvcount)+j*NZ+i, local[k*(NZ*recvcount)+j*NZ+i]);
			}
			printf("\n");
		}
		printf("\n\n");
	}
	
	if (taskId == 0)
	{
		MPI_Type_free(&type_ztPlane);
		MPI_Type_free(&type_ztPlane_resized);
	}
	MPI_Type_free(&type_ztxBlock);
	MPI_Type_free(&type_ztxBlock_resized);
	

	/*
	MPI_Type_vector(NT, NZ, NZ*NX, MPI_CHAR, &type_zfPlane);
	MPI_Type_commit(&type_zfPlane);
	int n = 1;
	local = (char*)calloc(n * NZ * NT, sizeof(char));
	if (taskId == 0)
		MPI_Send(global+2*NZ, n, type_zfPlane, 1, TAG, MPI_COMM_WORLD);
	else
		MPI_Recv(local, n * NZ * NT, MPI_CHAR, 0, TAG, MPI_COMM_WORLD, &status);
	if (taskId == 1)
	{
		printf("Printing local for rank %d\n", taskId);
		for (k = 0; k < NT; k++)
		{
			for (i = 0; i < NZ; i++)//j = 0; j < NX; j++
			{
				for (j = 0; j < n; j++)//i = 0; i < NZ; i++
				{
					printf("local[%d] = %c\t", k*(NZ*n)+j*NZ+i, local[k*(NZ*n)+j*NZ+i]);
				}
				printf("\n");
			}
			printf("\n\n");
		}
	}
	MPI_Type_free(&type_zfPlane);
	*/


    free(sendcounts);
    free(displs);
	free(local);
	MPI_Finalize();
	return 0;
}