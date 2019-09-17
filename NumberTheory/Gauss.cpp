const int MOD = 998244353;

int fpow(int a, int b) {
  int ans = 1;
  for(; b; b >>= 1, a = 1ll * a * a % MOD) {
    if(b & 1) ans = 1ll * ans * a % MOD;
  }
  return ans;
}

// make sure each element is in [0, MOD)
void Gauss(vector<vector<int>> &v) {
  const int m = v.size(), n = v[0].size();
  vector<int> id(n, -1);
  int i = 0, j = 0;
  while(i < m && j < n) {
    int mi = i;
    for(int k = i + 1; k < m; k++) {
      if(v[k][j] > v[mi][j]) {
        mi = k;
      }
    }
    if(v[mi][j] != 0) {
      if(i != mi) {
        for(int k = 0; k < n; k++) {
          swap(v[i][k], v[mi][k]);
        }
      }
      const int inv = fpow(v[i][j], MOD - 2);
      for(int k = j + 1; k < n; k++) {
        v[i][k] = 1ll * v[i][k] * inv % MOD;
      }
      id[j] = i;
      v[i][j] = 1;
      for(int r = i + 1; r < m; r++) {
        for(int c = j + 1; c < n; c++) {
          v[r][c] -= 1ll * v[r][j] * v[i][c] % MOD;
          if(v[r][c] < 0) v[r][c] += MOD;
        }
        v[r][j] = 0;
      }
      i++;
    }
    j++;
  }
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      if(v[i][j] && id[j] > i) {
        for(int k = j + 1; k < n; k++) {
          v[i][k] -= 1ll * v[i][j] * v[id[j]][k] % MOD;
          if(v[i][k] < 0) v[i][k] += MOD;
        }
        v[i][j] = 0;
      }
    }
  }
}
