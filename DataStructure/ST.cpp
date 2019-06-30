struct ST {
  vector<vector<int>> table;
  ST(vector<int> a = {}) {
    int n = a.size();
    table.resize(n, vector<int>(32 - __builtin_clz(n)));
    for (int i = 0; i < n; i++) {
      table[i][0] = a[i];
    }
    for (int j = 1; (1 << j) - 1 < n; j++) {
      for (int i = 0; i + (1 << j) - 1 < n; i++) {
        int x = table[i][j - 1], y = table[i + (1 << (j - 1))][j - 1];
        table[i][j] = min(x, y);
      }
    }
  }
  inline int getMin(int l, int r) {
    int k = 31 - __builtin_clz(r - l + 1);
    return min(table[l][k], table[r - (1 << k) + 1][k]);
  }
};

// ********************************************************************************************************
// Author: ConanYu
// faster than the above

struct ST {
  int *table, m, *lg;
  inline int idx(const int& i, const int& j) {
    return i * m + j;
  }
  ST(vector<int> a = {}) {
    int n = a.size();
    m = 32 - __builtin_clz(n);
    table = (int*) malloc(sizeof(int) * m * n);
    lg = (int*) malloc(sizeof(int) * (n + 1));
    lg[0] = lg[1] = 0;
    for(int i = 2; i <= n; i++) {
      lg[i] = lg[i >> 1] + 1;
    }
    for (int i = 0; i < n; i++) {
      table[idx(i, 0)] = a[i];
    }
    for (int j = 1; (1 << j) - 1 < n; j++) {
      for (int i = 0; i + (1 << j) - 1 < n; i++) {
        int x = table[idx(i, j - 1)], y = table[idx(i + (1 << (j - 1)), j - 1)];
        table[idx(i, j)] = min(x, y);
      }
    }
  }
  ~ST() {
    free(table);
    free(lg);
  }
  int getMin(int l, int r) {
    const int k = lg[r - l + 1];
    return min(table[idx(l, k)], table[idx(r - (1 << k) + 1, k)]);
  }
};
