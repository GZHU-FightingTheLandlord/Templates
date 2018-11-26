#include <bits/stdc++.h>
using namespace std;

const int maxn = 1e5 + 5;

vector<int> G[maxn];
int ver[maxn << 1], tot, First[maxn], dep[maxn << 1];
int dp[maxn << 1][20];
bool vis[maxn];

inline void init() {
    tot = 0;
    memset(vis, false, sizeof vis);
}

inline void addedge(int u, int v) {
    G[u].push_back(v);
    G[v].push_back(u);
}

void DFS(int u, int d) {
    vis[u] = true, ver[++tot] = u, First[u] = tot, dep[tot] = d;
    for (const int &v : G[u]) {
        if (!vis[v]) {
            DFS(v, d + 1);
            ver[++tot] = u, dep[tot] = d;
        }
    }
}

inline void ST() {
    for (int i = 1; i <= tot; i++) dp[i][0] = i;
    for (int j = 1; (1 << j) - 1<= tot; j++) {
        for (int i = 1; i + (1 << j) - 1 <= tot; i++) {
            int x = dp[i][j - 1], y = dp[i + (1 << (j - 1))][j - 1];
            dp[i][j] = (dep[x] < dep[y]) ? x : y;
        }
    }
}

inline int rmq(int l, int r) {
    int k = 31 - __builtin_clz(r - l + 1);
    // while ((1 << (k + 1)) <= r - l + 1) ++k;
    int x = dp[l][k], y = dp[r - (1 << k) + 1][k];
    return dep[x] < dep[y] ? x : y;
}

inline int lca(int u, int v) {
    int x = First[u], y = First[v];
    if (x > y) swap(x, y);
    return ver[rmq(x, y)];
}

// inline int lca(int u, int v) {
//     int l = First[u], r = First[v];
//     if (l > r) swap(l, r);
//     int k = 0;
//     while ((1 << (k + 1)) <= r - l + 1) ++k;
//     int x = dp[l][k], y = dp[r - (1 << k) + 1][k];
//     return ver[(dep[x] < dep[y]) ? x : y];
// }