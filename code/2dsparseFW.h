#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct cs_stable {
    int n;  // rows
    int m;  // cols
    int** atr;
} cs;

typedef struct dimx {
    int left;
    int right;
    int up;
    int down;
    int size0;
    int size1;
} dimensionx;

typedef struct etreeX {
    int h;
    int** Q;
    int* Qsize;
} etree_ds;

void assert(bool);
int Min(int, int);
void Swap(int*, int*);
void classicalFW(cs*);
void cs_malloc(cs*);
void cs_print_matrix(cs*, int);
bool find_array(int, int[], int);
int* matrix_to_1d(cs*, dimensionx);
void send_matrix(cs*, dimensionx, int, int, int, MPI_Comm);
void join_all_matrix(cs*, dimensionx, int, int, MPI_Comm);
int sum_a1(int, int);

void assert(bool err) {
    if (!err) {
        printf("Assert error!@!\n");
        exit(0);
    }
}

int Min(int x, int y) {
    if (x <= y) return x;
    return y;
}

void Swap(int* p, int* q) {
    int temp = *p;
    *p = *q;
    *q = temp;
}

void cs_malloc(cs* A) {
    A->atr = (int**)calloc(A->n, sizeof(int*));
    for (int i = 0; i < A->n; i++) {
        A->atr[i] = (int*)calloc(A->m, sizeof(int));
    }
}

void cs_print_matrix(cs* A, int my_cart_rank) {
    int n = A->n;
    int m = A->m;
    printf("P(%d) matrix\n", my_cart_rank);
    printf("---------------\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%d ", A->atr[i][j]);
        }
        printf("\n");
    }
}

void classicalFW(cs* A) {
    int n = A->n;
    assert(A->n == A->m);

    /*----- initializing -----*/
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
    //         A->atr[i][j] = A->atr[i][j];
    //         if (i == j) {
    //             A->atr[i][j] = 0;
    //         }
    //     }
    // }

    /*----- dynamic programming -----*/
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A->atr[i][j] = Min(A->atr[i][j], A->atr[i][k] + A->atr[k][j]);
            }
        }
    }
}

bool find_array(int key, int A[], int n) {
    bool flag = false;
    for (int i = 0; i < n; i++)
        if (A[i] == key) return true;
    return false;
}

int* matrix_to_1d(cs* A, dimensionx dimension) {
    int n = dimension.right - dimension.left + 1;
    int m = dimension.down - dimension.up + 1;
    int* p = (int*)calloc(n * m, sizeof(int));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            p[i * m + j] = A->atr[i + dimension.up][j + dimension.left];
        }
    }
    return p;
}

void send_matrix(cs* A, dimensionx dimension, int to_whom, int my_cart_rank,
                 int p, MPI_Comm comm) {
    int n = A->n;
    int m = A->m;
    int P[n * m];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            P[i * m + j] = A->atr[i][j];
        }
    }
    int left = dimension.left;
    int up = dimension.up;
    MPI_Send(&n, 1, MPI_INT, to_whom, 1, comm);
    MPI_Send(&m, 1, MPI_INT, to_whom, 2, comm);
    MPI_Send(&up, 1, MPI_INT, to_whom, 3, comm);
    MPI_Send(&left, 1, MPI_INT, to_whom, 4, comm);
    MPI_Send(P, n * m, MPI_INT, to_whom, 123, comm);
}

void join_all_matrix(cs* A, dimensionx dimension, int my_cart_rank, int p,
                     MPI_Comm comm) {
    MPI_Status status;
    int n, m, left, up;
    for (int i = 0; i < p; i++) {
        if (i == my_cart_rank) continue;
        MPI_Recv(&n, 1, MPI_INT, i, 1, comm, &status);
        MPI_Recv(&m, 1, MPI_INT, i, 2, comm, &status);
        int P[n * m];
        MPI_Recv(&up, 1, MPI_INT, i, 3, comm, &status);
        MPI_Recv(&left, 1, MPI_INT, i, 4, comm, &status);
        MPI_Recv(P, n * m, MPI_INT, i, 123, comm, &status);

        for (int u = 0; u < n; u++) {
            for (int v = 0; v < m; v++) {
                A->atr[u + up][v + left] = P[u * m + v];
            }
        }
    }
}

int sum_a1(int l, int r) {
    int sum = 0;
    for (int i = l; i <= r; i++) {
        sum += 1 << i;
    }
    return sum;
}
