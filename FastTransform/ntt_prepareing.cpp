// copy from https://codeforces.com/contest/1334/submission/76150604
namespace ntt {
  int base = 1, root = -1, max_base = -1;
  vector<int> rev = {0, 1}, roots = {0, 1};

  void init() {
    int temp = md - 1;
    max_base = 0;
    while (temp % 2 == 0) {
      temp >>= 1;
      ++max_base;
    }
    root = 2;
    while (true) {
      if (power(root, 1 << max_base) == 1 && power(root, 1 << max_base - 1) != 1) {
        break;
      }
      ++root;
    }
  }

  void ensure_base(int nbase) {
    if (max_base == -1) {
      init();
    }
    if (nbase <= base) {
      return;
    }
    assert(nbase <= max_base);
    rev.resize(1 << nbase);
    for (int i = 0; i < 1 << nbase; ++i) {
      rev[i] = rev[i >> 1] >> 1 | (i & 1) << nbase - 1;
    }
    roots.resize(1 << nbase);
    while (base < nbase) {
      int z = power(root, 1 << max_base - 1 - base);
      for (int i = 1 << base - 1; i < 1 << base; ++i) {
        roots[i << 1] = roots[i];
        roots[i << 1 | 1] = mul(roots[i], z);
      }
      ++base;
    }
  }

  void dft(vector<int> &a) {
    int n = a.size(), zeros = __builtin_ctz(n);
    ensure_base(zeros);
    int shift = base - zeros;
    for (int i = 0; i < n; ++i) {
      if (i < rev[i] >> shift) {
        swap(a[i], a[rev[i] >> shift]);
      }
    }
    for (int i = 1; i < n; i <<= 1) {
      for (int j = 0; j < n; j += i << 1) {
        for (int k = 0; k < i; ++k) {
          int x = a[j + k], y = mul(a[j + k + i], roots[i + k]);
          a[j + k] = (x + y) % md;
          a[j + k + i] = (x + md - y) % md;
        }
      }
    }
  }

  vector<int> multiply(vector<int> a, vector<int> b) {
    int need = a.size() + b.size() - 1, nbase = 0;
    while (1 << nbase < need) {
      ++nbase;
    }
    ensure_base(nbase);
    int sz = 1 << nbase;
    a.resize(sz);
    b.resize(sz);
    bool equal = a == b;
    dft(a);
    if (equal) {
      b = a;
    } else {
      dft(b);
    }
    int inv_sz = inv(sz);
    for (int i = 0; i < sz; ++i) {
      a[i] = mul(mul(a[i], b[i]), inv_sz);
    }
    reverse(a.begin() + 1, a.end());
    dft(a);
    a.resize(need);
    return a;
  }

  vector<int> inverse(vector<int> a) {
    int n = a.size(), m = n + 1 >> 1;
    if (n == 1) {
      return vector<int>(1, inv(a[0]));
    } else {
      vector<int> b = inverse(vector<int>(a.begin(), a.begin() + m));
      int need = n << 1, nbase = 0;
      while (1 << nbase < need) {
        ++nbase;
      }
      ensure_base(nbase);
      int sz = 1 << nbase;
      a.resize(sz);
      b.resize(sz);
      dft(a);
      dft(b);
      int inv_sz = inv(sz);
      for (int i = 0; i < sz; ++i) {
        a[i] = mul(mul(md + 2 - mul(a[i], b[i]), b[i]), inv_sz);
      }
      reverse(a.begin() + 1, a.end());
      dft(a);
      a.resize(n);
      return a;
    }
  }
}
