#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;

const int maxn = 1e4 + 5;

vector<int> e[maxn];
int dep[maxn], dp[maxn][15], maxb;

void init() {
    for (int i = 0; i < maxn; i++) {
        e[i].clear(), dep[i] = 0;
        memset(dp[i], -1, sizeof dp[i]);
    }
}

void DFS(int u, int d, int pre) {
    dp[u][0] = pre;
    dep[u] = d;
    for (int i = 0; i < (int)e[u].size(); i++) {
        if (e[u][i] != pre) {
            DFS(e[u][i], d + 1, u);
        }
    }
}

void run(int n) {
    maxb = 0;
    while ((1 << maxb) <= n) ++maxb;
    for (int j = 1; j < maxb; j++) {
        for (int i = 1; i <= n; i++) {
            (~dp[i][j - 1]) && (dp[i][j] = dp[dp[i][j - 1]][j - 1]);
        }
    }
}

int LCA(int u, int v) {
    if (dep[u] < dep[v]) swap(u, v);
    for (int j = maxb - 1; ~j; j--) {
        (dep[dp[u][j]] >= dep[v]) && (u = dp[u][j]);
    }
    if (u == v) return u;
    for (int j = maxb - 1; ~j; j--) {
        (dp[u][j] != dp[v][j]) && (u = dp[u][j], v = dp[v][j]);
    }
    return dp[u][0];
}