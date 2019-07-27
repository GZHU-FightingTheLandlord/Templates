struct LinearBasis {
  const static int MAXL = 50;
  long long a[MAXL + 1];
  LinearBasis() {
    memset(a, 0, sizeof a);
  }
  void insert(long long t) {
    for (int j = MAXL; j >= 0; j--) {
      if (!(t & (1ll << j))) {
        continue;
      }
      if (a[j]) {
        t ^= a[j];
      } else {
        for (int k = 0; k < j; k++) {
          if (t & (1ll << k)) {
            t ^= a[k];
          }
        }
        for (int k = j + 1; k <= MAXL; k++) {
          if (a[k] & (1ll << j)) {
            a[k] ^= t;
          }
        }
        a[j] = t;
        return;
      }
    }
  }
};


// 线性基区间最大值

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
