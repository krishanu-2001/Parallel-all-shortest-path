#include "2dsparseFW.h"

#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define INT_MAX (int)1e9

etree_ds data_etree;

int main(int argc, char *argv[]) {
    int ndims = 2, ierr;
    int p, my_rank, my_cart_rank;
    MPI_Comm comm2D;
    int dims[ndims], coord[ndims];
    int wrap_around[ndims];
    int reorder, nrows, ncols;
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Status status;
    MPI_Request req1;
    MPI_Request req2;
    int recv_count;
    double start, stop;

    /* process command line arguments*/
    if (argc == 3) {
        nrows = atoi(argv[1]);
        ncols = atoi(argv[2]);
        dims[0] = nrows; /* number of rows */
        dims[1] = ncols; /* number of columns */
        if ((nrows * ncols) != p) {
            if (my_rank == 0)
                // printf("ERROR: nrows*ncols)=%d * %d = %d != %d\n", nrows,
                // ncols,
                //        nrows * ncols, p);
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
        // printf("PW[%d]/[%d]: PEdims = [%d x %d] \n", my_rank, p, dims[0],
        //        dims[1]);
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
    // printf("PW[%d]: my_cart_rank PCM[%d], my coords = (%d,%d)\n", my_rank,
    //        my_cart_rank, coord[0], coord[1]);

    /*-- algorithm for elimination tree  --*/
    int h = log2(sqrt(p) + 1);
    data_etree.h = h;
    data_etree.Q = (int **)calloc(h, sizeof(int *));
    int elimination_levels = h;
    int level[(int)sqrt(p)];
    data_etree.Qsize = (int *)calloc(h, sizeof(int));
    for (int i = 0; i < h; i++) {
        int left = i * sqrt(p) / h;
        int right = (i + 1) * sqrt(p) / h - 1;
        if (i == h - 1) right = sqrt(p) - 1;
        data_etree.Qsize[i] = right - left + 1;
        data_etree.Q[i] = (int *)calloc(right - left + 1, sizeof(int));
        for (int j = left; j <= right; j++) {
            data_etree.Q[i][j - left] = j;
            level[j] = i;
        }
    }
    /*
        data_etree
        h = 2
        Q = [[0], [1]]
        Q_size = [1, 1]
    */

    /*----- input from file -----*/
    FILE *file;
    file = fopen("input", "r");
    int n;
    fscanf(file, "%d", &n);

    cs A;
    A.n = n;
    A.m = n;
    cs_malloc(&A);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%d", &A.atr[i][j]);
            if (A.atr[i][j] == 0 && i != j) A.atr[i][j] = INT_MAX;
        }
    }

    fclose(file);
    /*----- input completed -----*/

    /*
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
    */

    /*----- main function goes here -----*/
    // 1. eliminate supernodes
    /*
     _________         .    .
    (..       \_    ,  |\  /|
     \       O  \  /|  \ \/ /
      \______    \/ |   \  /
         vvvv\    \ |   /  |
         \^^^^  ==   \_/   |
          `\_   ===    \.  |
          / /\_   \ /      |
          |/   \_  \|      /
     snd         \________/
    */
    int r[2] = {n / dims[0], n / dims[1]};

    dimensionx dimension;
    int up = coord[0] * r[0];
    int down = (coord[0] == dims[0] - 1) ? n - 1 : (coord[0] + 1) * r[0] - 1;
    int left = coord[1] * r[1];
    int right = (coord[1] == dims[1] - 1) ? n - 1 : (coord[1] + 1) * r[1] - 1;

    dimension.up = up;
    dimension.left = left;
    dimension.right = right;
    dimension.down = down;
    dimension.size0 = (down - up + 1);
    dimension.size1 = (right - left + 1);

    cs iso;
    int k[2] = {down - up + 1, right - left + 1};
    iso.n = k[0];
    iso.m = k[1];
    cs_malloc(&iso);

    start = MPI_Wtime();

    // printf("\n\tLRUD - %d %d %d %d\n", left, right, up, down);

    for (int i = up; i <= down; i++) {
        for (int j = left; j <= right; j++) {
            iso.atr[i - up][j - left] = A.atr[i][j];
        }
    }

    int *Pkk;

    for (int l = 0; l < data_etree.h; l++) {
        // (i) diagonal update
        int Q_l_size = data_etree.Qsize[l];
        for (int q = 0; q < Q_l_size; q++) {
            if (coord[0] != data_etree.Q[l][q] ||
                coord[1] != data_etree.Q[l][q])
                continue;
            // 3. update R_l(1)

            classicalFW(&iso);

            for (int i = up; i <= down; i++) {
                for (int j = left; j <= right; j++) {
                    A.atr[i][j] = iso.atr[i - up][j - left];
                }
            }

            // (ii) panel update
            int n = dimension.size0;
            int m = dimension.size1;
            Pkk = (int *)calloc(n * m, sizeof(int));
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    Pkk[i * m + j] = iso.atr[i][j];
                }
            }

            for (int i = 0; i < dims[1]; i++) {
                if (i == coord[1]) continue;
                int rankv = dims[1] * coord[0] + i;
                MPI_Send(Pkk, dimension.size0 * dimension.size1, MPI_INT, rankv,
                         110, comm2D);
            }
            for (int i = 0; i < dims[1]; i++) {
                if (i == coord[0]) continue;
                int rankv = dims[1] * i + coord[0];
                MPI_Send(Pkk, dimension.size0 * dimension.size1, MPI_INT, rankv,
                         110, comm2D);
            }
        }

        /*----7. update R_l(2) ----*/
        /*
                ,-""""""-.
             /\j__/\  (  \`--.
        hjw  \`@_@'/  _)  >--.`.
            _{.:Y:_}_{{_,'    ) )
           {_}`-^{_} ```     (_/
        */

        for (int q = 0; q < Q_l_size; q++) {
            int k = data_etree.Q[l][q];
            if (coord[0] == coord[1]) continue;
            if (coord[0] == k) {
                // row
                int ranku = dims[1] * coord[0] + coord[0];
                Pkk = (int *)calloc(dimension.size0 * dimension.size1,
                                    sizeof(int));
                MPI_Recv(Pkk, dimension.size0 * dimension.size1, MPI_INT, ranku,
                         110, comm2D, &status);
                int n = dimension.size0;
                int m = dimension.size1;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        for (int t = 0; t < m; t++) {
                            iso.atr[i][j] = Min(iso.atr[i][j],
                                                Pkk[i * m + t] + iso.atr[t][j]);
                        }
                    }
                }
            } else if (coord[1] == k) {
                // col
                Pkk = (int *)calloc(dimension.size0 * dimension.size1,
                                    sizeof(int));
                int ranku = dims[1] * coord[1] + coord[1];
                MPI_Recv(Pkk, dimension.size0 * dimension.size1, MPI_INT, ranku,
                         110, comm2D, &status);
                int n = dimension.size0;
                int m = dimension.size1;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        for (int t = 0; t < m; t++) {
                            iso.atr[i][j] = Min(
                                iso.atr[i][j], +iso.atr[i][t] + Pkk[t * m + j]);
                        }
                    }
                }
            } else {
                continue;
            }
        }

        MPI_Barrier(comm2D);

        /*---- 9. R_l(2) BROADCAST to R_l(3) ----*/
        /*
        (             )
        `--(_   _)--'
            Y-Y
            /@@ \
            /     \
            `--'.  \             ,
                |   `.__________/)
        */
        int *Aik, *Akj;

        for (int q = 0; q < Q_l_size; q++) {
            int k = data_etree.Q[l][q];
            if (coord[0] == coord[1]) continue;
            if (coord[0] == k) {
                // row
                int n = dimension.size0;
                int m = dimension.size1;
                Aik = calloc(n * m, sizeof(int));
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        Aik[i * m + j] = iso.atr[i][j];
                    }
                }
                for (int i = 0; i < dims[0]; i++) {
                    if (i != coord[0]) {
                        int ranku = dims[1] * i + coord[1];
                        MPI_Send(Aik, dimension.size0 * dimension.size1,
                                 MPI_INT, ranku, 25, comm2D);
                    }
                }
            }
            if (coord[1] == k) {
                // col
                int n = dimension.size0;
                int m = dimension.size1;
                Akj = (int *)calloc(n * m, sizeof(int));
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        Akj[i * m + j] = iso.atr[i][j];
                    }
                }
                for (int i = 0; i < dims[1]; i++) {
                    if (i != coord[1]) {
                        int ranku = dims[1] * coord[0] + i;
                        MPI_Send(Akj, dimension.size0 * dimension.size1,
                                 MPI_INT, ranku, 35, comm2D);
                    }
                }
            }
        }

        /*---- 11. Update R_l(3) ----*/
        /*
                    __
                   /(`o
             ,-,  //  \\
            (,,,) ||   V
           (,,,,)\//
           (,,,/w)-'
           \,,/w)
           `V/uu
             / |
             | |
             o o
             \ |
        \,/  ,\|,.  \,/
        */

        for (int q = 0; q < Q_l_size; q++) {
            int k = data_etree.Q[l][q];
            if (coord[0] != k && coord[1] != k) {
                int ranku = dims[1] * k + coord[1];
                Aik = (int *)calloc(dimension.size0 * dimension.size1,
                                    sizeof(int));
                Akj = (int *)calloc(dimension.size0 * dimension.size1,
                                    sizeof(int));
                MPI_Recv(Akj, dimension.size0 * dimension.size1, MPI_INT, ranku,
                         25, comm2D, &status);
                ranku = dims[1] * coord[0] + k;
                MPI_Recv(Aik, dimension.size0 * dimension.size1, MPI_INT, ranku,
                         35, comm2D, &status);
                int n = dimension.size0;
                int m = dimension.size1;

                int left = dimension.left;
                int up = dimension.up;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        for (int t = 0; t < m; t++) {
                            iso.atr[i][j] = Min(
                                iso.atr[i][j], Aik[i * m + t] + Akj[t * m + j]);
                        }
                    }
                }
            }
        }

        /*---- 13. update R_l(4) - haaaaaaaaarrrrrrddddddd ----*/
        /*
                          ,,__
                ..  ..   / o._)                   .---.
               /--'/--\  \-'||        .----.    .'     '.
              /        \_/ / |      .'      '..'         '-.
            .'\  \__\  __.'.'     .'          i-._
              )\ |  )\ |      _.'
             // \\ // \\
            ||_  \\|_  \\_
        mrf '--' '--'' '--'
        */

        for (int q = 0; q < Q_l_size; q++) {
            int k = data_etree.Q[l][q];
            if (coord[0] == coord[1]) continue;
            if (coord[1] == k) {
                // Aik
                int n = dimension.size0;
                int m = dimension.size1;
                Aik = calloc(n * m, sizeof(int));
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        Aik[i * m + j] = iso.atr[i][j];
                    }
                }
                for (int a = l + 1; a < h; a++) {
                    int g = k - sum_a1(h - l + 1, h - 1);
                    for (int c = a; c < h; c++) {
                        int f = sum_a1(h + a - c, h - 1) + a - l;
                        assert(f < dims[0] && g < dims[1]);
                        int ranku = dims[1] * f + g;
                        if (my_cart_rank == ranku) continue;
                        // printf("\n\tP(%d) try to %d\n", my_cart_rank, ranku);
                        MPI_Send(Aik, dimension.size0 * dimension.size1,
                                 MPI_INT, ranku, 252, comm2D);
                        // printf("\n\tP(%d) sends %d\n", my_cart_rank, ranku);
                    }
                }
            }
            if (coord[0] == k) {
                // Akj
                int n = dimension.size0;
                int m = dimension.size1;
                Akj = (int *)calloc(n * m, sizeof(int));
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        Akj[i * m + j] = iso.atr[i][j];
                    }
                }
                for (int c = l + 1; c < h; c++) {
                    int g = k - sum_a1(h - l + 1, h - 1);
                    for (int a = l + 1; a < h; a++) {
                        int f = sum_a1(h + a - c, h - 1) + a - l;
                        assert(f < dims[0] && g < dims[1]);
                        int ranku = dims[1] * f + g;
                        if (my_cart_rank == ranku) continue;
                        // printf("\n\tP(%d) try to %d\n", my_cart_rank, ranku);
                        MPI_Send(Akj, dimension.size0 * dimension.size1,
                                 MPI_INT, ranku, 252, comm2D);
                        // printf("\n\tP(%d) sends to %d\n", my_cart_rank,
                        // ranku);
                    }
                }
            }
        }

        /*---- 18. pfk to R_l(4) ----*/
        /*
        ░░░░░░░░░░░░░░░░░░░░░▄▀░░▌
        ░░░░░░░░░░░░░░░░░░░▄▀▐░░░▌
        ░░░░░░░░░░░░░░░░▄▀▀▒▐▒░░░▌
        ░░░░░▄▀▀▄░░░▄▄▀▀▒▒▒▒▌▒▒░░▌
        ░░░░▐▒░░░▀▄▀▒▒▒▒▒▒▒▒▒▒▒▒▒█
        ░░░░▌▒░░░░▒▀▄▒▒▒▒▒▒▒▒▒▒▒▒▒▀▄
        ░░░░▐▒░░░░░▒▒▒▒▒▒▒▒▒▌▒▐▒▒▒▒▒▀▄
        ░░░░▌▀▄░░▒▒▒▒▒▒▒▒▐▒▒▒▌▒▌▒▄▄▒▒▐
        ░░░▌▌▒▒▀▒▒▒▒▒▒▒▒▒▒▐▒▒▒▒▒█▄█▌▒▒▌
        ░▄▀▒▐▒▒▒▒▒▒▒▒▒▒▒▄▀█▌▒▒▒▒▒▀▀▒▒▐░░░▄
        ▀▒▒▒▒▌▒▒▒▒▒▒▒▄▒▐███▌▄▒▒▒▒▒▒▒▄▀▀▀▀
        ▒▒▒▒▒▐▒▒▒▒▒▄▀▒▒▒▀▀▀▒▒▒▒▄█▀░░▒▌▀▀▄▄
        ▒▒▒▒▒▒█▒▄▄▀▒▒▒▒▒▒▒▒▒▒▒░░▐▒▀▄▀▄░░░░▀
        ▒▒▒▒▒▒▒█▒▒▒▒▒▒▒▒▒▄▒▒▒▒▄▀▒▒▒▌░░▀▄
        ▒▒▒▒▒▒▒▒▀▄▒▒▒▒▒▒▒▒▀▀▀▀▒▒▒▄▀
        */
        int n = dimension.size0;
        int m = dimension.size1;
        int TOT[n * m];
        int P[n * m];
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < m; v++) {
                TOT[u * m + v] = iso.atr[u][v];
            }
        }
        for (int q = 0; q < Q_l_size; q++) {
            int k = data_etree.Q[l][q];
            for (int i = k + 1; i < h; i++) {
                for (int j = i; j < h; j++) {
                    // a = level[i];
                    // c = level[j];
                    int a = level[i];
                    int c = level[j];
                    // 23. reduce
                    for (int q1 = 0; q1 < Q_l_size; q1++) {
                        int k1 = data_etree.Q[l][q1];
                        if (k1 < i) {
                            // D(i)
                            int f = sum_a1(h + a - c, h - 1) + a - l;
                            int g = k - sum_a1(h - l + 1, h - 1);
                            assert(f < dims[0] && g < dims[1]);
                            if (coord[0] == f && coord[1] == g) {
                                // printf("\n\tkfg %d %d %d\n", k, f, g);
                                // printf("\n\tk1ij %d %d %d\n", k1, i, j);
                                int ranku = dims[1] * i + k1;
                                if (my_cart_rank != ranku && i == j) {
                                    // printf("\n\tP(%d) try recv from %d\n",
                                    //        my_cart_rank, ranku);
                                    MPI_Recv(Aik, n * m, MPI_INT, ranku, 252,
                                             comm2D, &status);
                                    // printf("\n\tP(%d) recv %d\n",
                                    // my_cart_rank,
                                    //        ranku);
                                }

                                int rankv = dims[1] * k1 + j;
                                if (my_cart_rank != rankv && i == j) {
                                    // printf("\n\tP(%d) try recv to %d\n",
                                    //        my_cart_rank, rankv);
                                    MPI_Recv(Akj, n * m, MPI_INT, rankv, 252,
                                             comm2D, &status);
                                    // printf("\n\tP(%d) recv %d\n",
                                    // my_cart_rank,
                                    //        rankv);
                                }
                                for (int u = 0; u < n; u++) {
                                    for (int v = 0; v < m; v++) {
                                        P[u * m + v] = INT_MAX;
                                        for (int t = 0; t < m; t++) {
                                            P[u * m + v] =
                                                Min(P[u * m + v],
                                                    Aik[u * m + t] +
                                                        Akj[t * m + v]);
                                        }
                                    }
                                }
                                // 23. pfg reduces to pij
                                int rankij = dims[1] * i + j;
                                if (my_cart_rank != rankij) {
                                    MPI_Send(P, n * m, MPI_INT, rankij, 119,
                                             comm2D);
                                }
                            }
                            if (coord[0] == i && coord[1] == j) {
                                // 23. pfg reduces to pij
                                int ranku = dims[1] * f + g;
                                if (my_cart_rank != ranku) {
                                    MPI_Recv(P, n * m, MPI_INT, ranku, 119,
                                             comm2D, &status);
                                }
                                for (int y1 = 0; y1 < n * m; y1++) {
                                    TOT[y1] = Min(TOT[y1], P[y1]);
                                }
                            }
                        }
                    }

                    // 25. symmetic matrix
                    if (coord[0] == i && coord[1] == j) {
                        for (int u = 0; u < n; u++) {
                            for (int v = 0; v < m; v++) {
                                iso.atr[u][v] = TOT[u * m + v];
                            }
                        }
                        int rankji = dims[1] * j + i;
                        if (i != j)
                            MPI_Send(TOT, n * m, MPI_INT, rankji, 432, comm2D);
                    }

                    if (coord[0] == j && coord[1] == i && i != j) {
                        int rankij = dims[1] * i + j;
                        MPI_Recv(TOT, n * m, MPI_INT, rankij, 432, comm2D,
                                 &status);
                        for (int u = 0; u < n; u++) {
                            for (int v = 0; v < m; v++) {
                                iso.atr[u][v] = TOT[u * m + v];
                            }
                        }
                    }
                }
            }
        }
        // the end
    }
    /*
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done
      It always seems impossible  until it's done

    */
    MPI_Barrier(comm2D);
    stop = MPI_Wtime();
    if (my_rank != 0) {
        send_matrix(&iso, dimension, 0, my_cart_rank, p, comm2D);
    }

    if (my_rank == 0) {
        printf("\nRuntime : %f (us)\n\n", (stop - start) * 1e6);

        file = fopen("output_sparse", "w");

        int n = dimension.size0;
        int m = dimension.size1;
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < m; v++) {
                A.atr[u + dimension.up][v + dimension.left] = iso.atr[u][v];
            }
        }

        join_all_matrix(&A, dimension, my_cart_rank, p, comm2D);

        /*----- print shortest path matrix ----*/
        n = A.n;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fprintf(file, "%12d ", A.atr[i][j]);
            }
            fprintf(file, "\n");
        }
        printf("Sparse Floyd-Warshall completed!\n");
    }
    MPI_Comm_free(&comm2D);
    MPI_Finalize();
    return 0;
}