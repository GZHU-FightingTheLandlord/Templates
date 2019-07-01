namespace cipolla {
  long long p, w;
  struct cpx {
    long long x, y;
    cpx(long long xx = 0, long long yy = 0) {
      x = xx, y = yy;
    }
    friend cpx operator * (const cpx & a, const cpx & b) {
      return cpx((a.x * b.x + a.y * b.y % p * w % p) % p, (a.x * b.y + a.y * b.x) % p);
    }
    friend cpx operator % (const cpx & a, const long long & b) {
      return a;
    }
  };
  template<typename T>
  T qpow(T a, long long b) {
    T ans(1);
    for(; b; a = a * a % p, b >>= 1) {
      if(b & 1) {
        ans = 1ll * ans * a % p;
      }
    }
    return ans;
  }
  // p >= 3 if return `a`, then `p - a` is the ans too.
  long long solve(long long n, long long pp) {
    p = pp, n = ((n % p) + p) % p;
    long long q = qpow(n, (p - 1) >> 1);
    if(q == 0 || q == p - 1) {
      return -1;
    }
    for(q = rand() % p; qpow((q * q - n + p) % p, (p - 1) >> 1) != p - 1; q = rand() % p);
    cpx t(q, 1);
    w = (q * q - n + p) % p;
    return qpow(t, (p + 1) >> 1).x;
  }
}
