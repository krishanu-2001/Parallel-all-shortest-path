/* mpi-cart-2D.c -- test basic -cartesian functions
 * Written by Mary Thomas- Updated Mar, 2015
 * Based loosely on code from Pacheco’97,
 * Chap 7, pp. 121 & ff in PPMPI */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int ndims = 2, ierr;
    int p, my_rank, my_cart_rank;
    MPI_Comm comm2D;
    int dims[ndims], coord[ndims];
    int wrap_around[ndims];
    int reorder, nrows, ncols;
    /* start up initial MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* process command line arguments*/
    if (argc == 3) {
        nrows = atoi(argv[1]);
        ncols = atoi(argv[2]);
        dims[0] = nrows; /* number of rows */
        dims[1] = ncols; /* number of columns */
        if ((nrows * ncols) != p) {
            if (my_rank == 0)
                printf("ERROR: nrows*ncols)=%d * %d = %d != %d\n", nrows, ncols,
                       nrows * ncols, p);
            MPI_Finalize();
            exit(0);
        }
    } else {
        nrows = ncols = (int)sqrt(p);
        dims[0] = dims[1] = 0;
    }
    /* create cartesian topology for processes */
    MPI_Dims_create(p, ndims, dims);
    if (my_rank == 0)
        printf("PW[%d]/[%d]: PEdims = [%d x %d] \n", my_rank, p, dims[0],
               dims[1]);
    /* create cartesian mapping */
    wrap_around[0] = wrap_around[1] = 0;  // set periodicity
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder,
                           &comm2D);
    if (ierr != 0) printf("ERROR[%d] creating CART\n", ierr);
    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm2D, my_rank, ndims, coord);
    /* use my coords to find my rank in cartesian group*/
    MPI_Cart_rank(comm2D, coord, &my_cart_rank);
    printf("PW[%d]: my_cart_rank PCM[%d], my coords = (%d,%d)\n", my_rank,
           my_cart_rank, coord[0], coord[1]);
    MPI_Comm_free(&comm2D);
    MPI_Finalize();
} /* main */