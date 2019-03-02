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