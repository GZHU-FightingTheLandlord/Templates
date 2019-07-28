typedef unsigned long long ull;
struct LB {
  const static int L = 64; // insert(x) 0 <= x < (1 << L)
  ull a[L];
  LB() { this->init(); }
  void init() { memset(a, 0, sizeof a); }
  ull &operator[](const size_t &id) { return a[id]; }
  const ull &operator[](const size_t &id) const { return a[id]; }
  bool insert(ull x) {
    for(int i = L - 1; ~i; i--) {
      if((x >> i) & 1) {
        if(!a[i]) {
          a[i] = x;
          return true;
        } else {
          x ^= a[i];
        }
      }
      if(!x) {
        return false;
      }
    }
    return true;
  }
  // 高斯消元 化为行简化阶梯型
  void operator()() {
    for(int i = L - 1; ~i; i--) {
      if(a[i]) {
        for(int j = i + 1; j < L; j++) {
          a[j] = min(a[j], a[j] ^ a[i]);
        }
      }
    }
  }
  // 交
  friend LB operator&(const LB &A, const LB &B) {
    LB C, D, E;
    for(int i = L - 1; ~i; i--) {
      if(A[i]) {
        C.insert(A[i]);
        D[i] = 1ull << i;
      }
    }
    for(int i = 0; i < L; i++) {
      if(!B[i]) {
        continue;
      }
      bool can = true;
      ull v = 0, x = B[i];
      for(int j = L - 1; ~j; j--) {
        if((x >> j) & 1) {
          if(C[j]) {
            x ^= C[j], v ^= D[j];
          } else {
            can = false, C[j] = x, D[j] = v;
            break;
          }
        }
      }
      if(can) {
        ll m = 0;
        for(int j = L - 1; ~j; j--) {
          if((v >> j) & 1) {
            m ^= A[j];
          }
        }
        E.insert(m);
      }
    }
    return E;
  }
  // 并
  friend LB operator|(const LB &x, const LB &y) {
    LB z;
    for(int i = 0; i < L; i++) if(x[i]) z.insert(x[i]);
    for(int i = 0; i < L; i++) if(y[i]) z.insert(y[i]);
    return z;
  }
};