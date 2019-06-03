namespace fwt {

  inline void fwt(int &a, int &b) {
    int x = a, y = b;
    a = x + y, b = x - y;
    // xor: a = x + y, b = x - y
    // and: a = x + y
    // or: b = x + y
    // 记得取模
  }

  inline void ifwt(int &a, int &b) {
    int x = a, y = b;
    b = y - x;
    a = (x + y) >> 1, b = (x - y) >> 1;
    // xor: a = (x + y) / 2, b = (x - y) / 2
    // and: a = x - y
    // or: a = y - x
    // 记得取模
  }

  void trans(vector<int> &a, int n, bool is) {
    for (int d = 1; d < n; d <<= 1) {
      for (int i = 0, k = d << 1; i < n; i += k) {
        for (int j = 0; j < d; j++) {
          if(is) fwt(a[i + j], a[i + j + d]);
          else ifwt(a[i + j], a[i + j + d]);
        }
      }
    }
  }
  
  // attendtion: 该函数会修改a和b且答案存储在a中
  void solve(vector<int> &a, vector<int> &b) { 
    int n = max(a.size(), b.size());
    n = 1 << (32 - __builtin_clz(n) - (n - (n & (-n)) == 0));
    a.resize(n);
    b.resize(n);
    trans(a, n, true);
    trans(b, n, true);
    for(int i = 0; i < n; i++) {
      a[i] = a[i] * b[i]; // 记得取模
    }
    trans(a, n, false);
  }
}
