/* File:      mpi_floyd.c
 * Purpose:   Implement a parallel version of Floyd's algorithm for finding
 *            the least cost path between each pair of vertices in a labelled
 *            digraph.
 *
 * Compile:   mpicc -g -Wall -o mpi_floyd mpi_floyd.c
 * Run:       mpiexec -n <number of processes> ./mpi_floyd
 *
 * Input:     n, the number of vertices
 *            mat, the adjacency matrix
 * Output:    mat, after being updated by floyd so that it contains the
 *            costs of the cheapest paths between all pairs of vertices.
 *
 * Notes:
 * 1.  n, the number of vertices should be evenly divisible by p, the
 *     number of processes.
 * 2.  The entries in the matrix should be nonnegative ints:  0 on the
 *     diagonal, positive off the diagonal.  Infinity should be indicated
 *     by the constant INFINITY.  (See below.)
 * 3.  The matrix is distributed by block rows.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for debugging */

const int INFINITY = (int)1e9;

void Read_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm);
void Print_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm);
void Floyd(int local_mat[], int n, int my_rank, int p, MPI_Comm comm);
int Owner(int k, int p, int n);
void Copy_row(int local_mat[], int n, int p, int row_k[], int k);
void Print_row(int local_mat[], int n, int my_rank, int i);

int main(int argc, char* argv[]) {
    int n;
    int* local_mat;
    MPI_Comm comm;
    int p, my_rank;

    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    if (my_rank == 0) {
        FILE* file;
        file = fopen("input", "r");
        int n;
        fscanf(file, "%d", &n);
        fclose(file);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    local_mat = malloc(n * n / p * sizeof(int));

    if (my_rank == 0) printf("Enter the local_matrix\n");
    Read_matrix(local_mat, n, my_rank, p, comm);
    // if (my_rank == 0) printf("We got\n");
    // Print_matrix(local_mat, n, my_rank, p, comm);
    // if (my_rank == 0) printf("\n");

    double start, end;

    MPI_Barrier(comm); /* IMPORTANT */
    start = MPI_Wtime();

    Floyd(local_mat, n, my_rank, p, comm);

    MPI_Barrier(comm); /* IMPORTANT */
    end = MPI_Wtime();

    if(my_rank == 0){
        printf("Runtime : %f (us)\n\n", (end - start) * 1e6);
    }

    // if (my_rank == 0) printf("The solution is:\n");
    // Print_matrix(local_mat, n, my_rank, p, comm);

    free(local_mat);
    MPI_Finalize();

    return 0;
} /* main */

/*---------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in the local_matrix on process 0 and scatter it using a
 *            block row distribution among the processes.
 * In args:   All except local_mat
 * Out arg:   local_mat
 */
void Read_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm) {
    int i, j;
    int* temp_mat = NULL;

    if (my_rank == 0) {
        FILE* file;
        file = fopen("input", "r");
        int n_extra;
        fscanf(file, "%d", &n_extra);
        temp_mat = malloc(n * n * sizeof(int));
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                int x;
                fscanf(file, "%d", &x);
                if (x == 0 && i != j) x = INFINITY;
                temp_mat[i * n + j] = x;
            }
        fclose(file);

        MPI_Scatter(temp_mat, n * n / p, MPI_INT, local_mat, n * n / p, MPI_INT,
                    0, comm);
        free(temp_mat);
    } else {
        MPI_Scatter(temp_mat, n * n / p, MPI_INT, local_mat, n * n / p, MPI_INT,
                    0, comm);
    }

} /* Read_matrix */

/*---------------------------------------------------------------------
 * Function:  Print_row
 * Purpose:   Convert a row of the matrix to a string and then print
 *            the string.  Primarily for debugging:  the single string
 *            is less likely to be corrupted when multiple processes
 *            are attempting to print.
 * In args:   All
 */
void Print_row(int local_mat[], int n, int my_rank, int i) {
    char char_int[100];
    char char_row[1000];
    int j, offset = 0;

    for (j = 0; j < n; j++) {
        if (local_mat[i * n + j] == INFINITY)
            sprintf(char_int, "i ");
        else
            sprintf(char_int, "%d ", local_mat[i * n + j]);
        sprintf(char_row + offset, "%s", char_int);
        offset += strlen(char_int);
    }
    printf("Proc %d > row %d = %s\n", my_rank, i, char_row);
} /* Print_row */

/*---------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Gather the distributed matrix onto process 0 and print it.
 * In args:   All
 */
void Print_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm) {
    int i, j;
    int* temp_mat = NULL;

    if (my_rank == 0) {
        temp_mat = malloc(n * n * sizeof(int));
        MPI_Gather(local_mat, n * n / p, MPI_INT, temp_mat, n * n / p, MPI_INT,
                   0, comm);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                if (temp_mat[i * n + j] == INFINITY)
                    printf("i ");
                else
                    printf("%d ", temp_mat[i * n + j]);
            printf("\n");
        }
        free(temp_mat);
    } else {
        MPI_Gather(local_mat, n * n / p, MPI_INT, temp_mat, n * n / p, MPI_INT,
                   0, comm);
    }
} /* Print_matrix */

/*---------------------------------------------------------------------
 * Function:    Floyd
 * Purpose:     Implement a distributed version of Floyd's algorithm for
 *              finding the shortest path between all pairs of vertices.
 *              The adjacency matrix is distributed by block rows.
 * In args:     All except local_mat
 * In/out arg:  local_mat:  on input the adjacency matrix.  On output
 *              the matrix of lowests costs between all pairs of
 *              vertices
 */
void Floyd(int local_mat[], int n, int my_rank, int p, MPI_Comm comm) {
    int global_k, local_i, global_j, temp;
    int root;
    int* row_k = malloc(n * sizeof(int));

    for (global_k = 0; global_k < n; global_k++) {
        root = Owner(global_k, p, n);
        if (my_rank == root) Copy_row(local_mat, n, p, row_k, global_k);
        MPI_Bcast(row_k, n, MPI_INT, root, comm);
        for (local_i = 0; local_i < n / p; local_i++)
            for (global_j = 0; global_j < n; global_j++) {
                temp = local_mat[local_i * n + global_k] + row_k[global_j];
                if (temp < local_mat[local_i * n + global_j])
                    local_mat[local_i * n + global_j] = temp;
            }
    }
    free(row_k);
} /* Floyd */

/*---------------------------------------------------------------------
 * Function:  Owner
 * Purpose:   Return rank of process that owns global row k
 * In args:   All
 */
int Owner(int k, int p, int n) { return k / (n / p); } /* Owner */

/*---------------------------------------------------------------------
 * Function:  Copy_row
 * Purpose:   Copy the row with *global* subscript k into row_k
 * In args:   All except row_k
 * Out arg:   row_k
 */
void Copy_row(int local_mat[], int n, int p, int row_k[], int k) {
    int j;
    int local_k = k % (n / p);

    for (j = 0; j < n; j++) row_k[j] = local_mat[local_k * n + j];
} /* Copy_row */