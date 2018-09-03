#include <bits/stdc++.h>
using namespace std;

const int MAX = 1e6 + 5;

int dp[MAX][25];

void rmq(int n) {
    int len = (int)(log(n) / log(2.0));
    for (int j = 1; j <= len; j++) {
        for (int i = 1; i + (1 << j) - 1 <= n; i++) {
            dp[i][j] = min(dp[i][j - 1], dp[i + (1 << (j - 1))][j - 1]);
        }
    }
}

int query(int l, int r) {
    int p = (int)(log(r - l + 1) / log(2.0));
    return min(dp[l][p], dp[r - (1 << p) + 1][p]);
}
