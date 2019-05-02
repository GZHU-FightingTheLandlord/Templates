// get this code from "https://codeforces.com/contest/1096/submission/47636051"
// author: palayutm

namespace ntt {
  int qpow(int a, int t, int mod) {
    ll b = 1;
    for (; t; t >>= 1, a = (ll)a * a % mod) {
      if (t & 1) b = b * a % mod; 
    }
    return b;
  }
  int revv(int x, int bits) {
    int ret = 0;
    for (int i = 0; i < bits; i++) {
      ret <<= 1, ret |= x & 1, x >>= 1;
    }
    return ret;
  }
  void ntt(vector<int> &a, bool rev, int mod, int root) {
    int n = (int)a.size(), bits = 31 - __builtin_clz(n);
    for (int i = 0; i < n; i++) {
      int j = revv(i, bits);
      if (i < j) swap(a[i], a[j]);
    }
    for (int k = 1; k < n; k <<= 1) {
      int e = qpow(root, (mod - 1) / 2 / k, mod);
      if (rev) e = qpow(e, mod - 2, mod);
      for (int i = 0; i < n; i += 2 * k) {
        ll w = 1;
        for (int j = 0; j < k; j++, w = w * e % mod) {
          int x = a[i + j], y = w * a[i + j + k] % mod;
          a[i + j] = (x + y) % mod, a[i + j + k] = (x - y + mod) % mod;
        }
      }
    }
    if (rev) {
      int inv = qpow(n, mod - 2, mod);
      for (int i = 0; i < n; i++) a[i] = 1ll * a[i] * inv % mod;
    }
  }
  // mod = 998244353 = (119 << 23) + 1, root = 3, // = (119 << 23, 3)
  // For p < 2^30, (5 << 25, 3), (7 << 26, 3),
  // (479 << 21, 3) and (483 << 21, 5), last two are > 10^9
  vector<int> conv(const vector<int>& a, const vector<int>& b,\
                  const int mod = (119 << 23) + 1, int root = 3) {
    int sz = (int)a.size() + (int)b.size() - 1;
    int L = sz > 1 ? (32 - __builtin_clz(sz - 1)) : 0, n = 1 << L;
    vector<int> av(n), bv(n);
    copy(a.begin(), a.end(), av.begin());
    copy(b.begin(), b.end(), bv.begin());
    ntt(av, false, mod, root), ntt(bv, false, mod, root);
    for (int i = 0; i < n; i++) {
      av[i] = 1ll * av[i] * bv[i] % mod;
    }
    ntt(av, true, mod, root);
    av.resize(sz);
    return av;
  }
}
