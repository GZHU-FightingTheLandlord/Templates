const int N = 505, M = 10;
int ST[N][N][M][M], a[N][N], lg[N];

void build(int n, int m) {
  for(int i = 2; i <= max(n, m); i++) {
    lg[i] = lg[i / 2] + 1;
  }
  for(int i = 1; i <= n; i++) {
    for(int j = 1; j <= m; j++) {
      ST[i][j][0][0] = a[i][j];
    }
  }
  for(int i = 1; i <= n; i++) {
    for(int k = 1; (1 << k) <= m; k++) {
      for(int j = 1; j + (1 << k) - 1 <= m; j++) {
        ST[i][j][0][k] = max(ST[i][j][0][k - 1], ST[i][j + (1 << (k - 1))][0][k - 1]);
      }
    }
  }
  for(int k1 = 1; (1 << k1) <= n; k1++) {
    for(int i = 1; i + (1 << k1) - 1 <= n; i++) {
      for(int k2 = 0; (1 << k2) <= m; k2++) {
        for(int j = 1; j + (1 << k2) - 1 <= m; j++) {
          ST[i][j][k1][k2] = max(ST[i][j][k1 - 1][k2], ST[i + (1 << (k1 - 1))][j][k1 - 1][k2]);
        }
      }
    }
  }
}

int query(int x1, int y1, int x2, int y2) {
  if(x2 < x1 || y2 < y1) return 0;
  const int k1 = lg[x2 - x1 + 1], k2 = lg[y2 - y1 + 1];
  x2 = x2 - (1 << k1) + 1;
  y2 = y2 - (1 << k2) + 1;
  return max({ST[x1][y1][k1][k2], ST[x1][y2][k1][k2], ST[x2][y1][k1][k2], ST[x2][y2][k1][k2]});
}
