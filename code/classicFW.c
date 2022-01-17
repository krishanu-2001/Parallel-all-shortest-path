#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define INT_MAX (int)1e9

/*
 * classicalFW ASPS
 * Time - O(n^3)
 * Space - O(n^2)
 */

int Min(int x, int y) {
    if (x <= y) return x;
    return y;
}

int main() {
    /*----- input from file -----*/
    FILE *file;
    file = fopen("input", "r");
    int n;
    fscanf(file, "%d", &n);
    int adj[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(file, "%d", &adj[i][j]);
            if (adj[i][j] == 0) adj[i][j] = INT_MAX;
        }
    }
    fclose(file);
    /*----- input completed -----*/

    clock_t t;
    t = clock();

    int dp[n][n];

    /*----- initializing -----*/
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dp[i][j] = adj[i][j];
            if (i == j) {
                dp[i][j] = 0;
            }
        }
    }

    /*----- dynamic programming -----*/
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dp[i][j] = Min(dp[i][j], dp[i][k] + dp[k][j]);
            }
        }
    }

    t = clock() - t;
    double time_taken = ((double)t) / CLOCKS_PER_SEC;  // in seconds
    printf("\n ClassicalFW Runtime took %f (us)\n", time_taken * 1e6);

    file = fopen("output_classic", "w");
    /*----- print shortest path matrix ----*/
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%12d ", dp[i][j]);
        }
        fprintf(file, "\n");
    }
    printf("\nClassical Floyd-Warshall completed!\n");

    return 0;
}