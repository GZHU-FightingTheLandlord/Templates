typedef unsigned long long ull;
struct LB {
  const static int L = 64; // insert(x) 0 <= x < (1 << L)
  ull a[L];
  LB() { this->init(); }
  void init() { memset(a, 0, sizeof a); }
  ull &operator[](const size_t &id) { return a[id]; }
  const ull &operator[](const size_t &id) const { return a[id]; }
  // 询问x是否在线性基中可以仿造insert函数来写
  // 即将`return true;`上面三行删去 然后把返回值取反
  bool insert(ull x) {
    for(int i = L - 1; ~i; i--) {
      if((x >> i) & 1) {
        if(!a[i]) {
          for(int j = 0; j < i; j++) if((x >> j) & 1) x ^= a[j];
          for(int j = i + 1; j < L; j++) if((a[j] >> i) & 1) a[j] ^= x;
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

// 线性基区间最大值
// hdu-6579 https://vjudge.net/problem/HDU-6579
// 还需改进

namespace LBRMQ {
  const int N = 1e6 + 10, L = 32;
  int b[N][L], pre[N][L];
  void init() {	
    memset(b[0], 0, sizeof b[0]);	
    memset(pre[0], 0, sizeof pre[0]);	
  }	
  // index start from 1	
  void add(int x, int r) {	
    int maxl = r;	
    memcpy(b[r], b[r - 1], sizeof(int) * L);	
    memcpy(pre[r], pre[r - 1], sizeof(int) * L);	
    for(int i = L - 1; ~i; i--) {	
      if((x >> i) & 1) {	
        if(!b[r][i]) {	
          b[r][i] = x;	
          pre[r][i] = maxl;	
          return;	
        }	
        if(pre[r][i] < maxl) {	
          swap(pre[r][i], maxl);	
          swap(b[r][i], x);	
        }	
        x ^= b[r][i];	
      }	
    }	
  }	
  int query(int l, int r) {	
    int ans = 0;	
    for(int i = L - 1; ~i; i--) {	
      if(pre[r][i] >= l) {	
        ans = max(ans, ans ^ b[r][i]);	
      }	
    }	
    return ans;	
  }
}
