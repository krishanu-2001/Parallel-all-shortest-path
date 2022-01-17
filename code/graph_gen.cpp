#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

#define pb push_back
#define ff first
#define ss second

void fill_adjlist(vector<pair<ll, ll>> adj[], ll n) {
    for (int i = 1; i <= n; i++) {
        for (int j = 0; j <= 32; j++) {
            if (((1LL << j) & i)) {
                ll w = i + (1LL << j);
                adj[i].pb({(1LL << j), w});
            }
        }
    }
}

int main(int argc, char** argv) {
#ifndef ONLINE_JUDGE
    freopen("input", "w", stdout);
#endif  // ONLINE_JUDGE

    int n = atoi(argv[1]);
    int mat[n + 1][n + 1];
    memset(mat, 0, sizeof(mat));
    vector<pair<ll, ll>> adj[n + 1];

    fill_adjlist(adj, n);

    for (int i = 1; i <= n; i++) {
        for (auto e : adj[i]) {
            int u = i;
            int v = e.ff;
            int w = e.ss;
            if (u == v) w = 0;
            mat[u][v] = w;
            mat[v][u] = w;
        }
    }

    cout << n << "\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            cout << mat[i][j] << " ";
        }
        cout << "\n";
    }

    return 0;
}
