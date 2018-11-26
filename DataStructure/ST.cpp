const int maxn = 1e5 + 5;
const int maxm = 20; // larger than log2(maxn)

int n;
int arr[maxn], dp[maxn][maxm];

void ST() {
    for (int i = 1; i <= n; i++) dp[i][0] = arr[i];
    for (int j = 1; (1 << j) - 1 <= n; j++) {
        for (int i = 1; i + (1 << j) - 1 <= n; i++) {
            dp[i][j] = max(dp[i][j - 1], dp[i + (1 << (j - 1))][j - 1]);
        }
    }
}

int query(int l, int r) {
    int k = 31 - __builtin_clz(r - l + 1);
    // while ((1 << (k + 1)) <= r - l + 1) ++k;
    return max(dp[l][k], dp[r - (1 << k) + 1][k]);
}