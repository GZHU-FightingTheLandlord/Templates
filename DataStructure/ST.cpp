const int N = 1e5 + 5, M = 18;

int dp[N][M], lg[N];

void build(int *arr, int n) {
  lg[0] = -1;
  for (int i = 1; i <= n; i++) {
    dp[i][0] = arr[i], lg[i] = lg[i >> 1] + 1;
  }
  for (int j = 1; (1 << j) <= n; j++) {
    for (int i = 1; i + (1 << j) - 1 <= n; i++) {
      int x = dp[i][j - 1], y = dp[i + (1 << (j - 1))][j - 1];
      dp[i][j] = min(x, y);
    }
  }
}

int getMin(int l, int r) {
  int k = lg[r - l + 1];
  return min(dp[l][k], dp[r - (1 << k) + 1][k]);
}