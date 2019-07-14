// 引用 exgcd 和 fft::multiply_mod
// 确保所有输入的数在[0, MOD)的区间中
// MOD < 1073741823 以及 最好是质数
// 传入的vector为{a_0, a_1, a_2, ..., a_n} 即认定为 y=\sum_{i=0}^{n}a_i\cdot x^i
namespace poly {
  const int MOD = 998244353ll;
  vector<int> inv(const vector<int> &a) {
    if(a.size() == 1) {
      const int inv = exgcd(a[0], MOD).first;
      return vector<int>(1, inv < 0 ? inv + MOD : inv);
    }
    const int na = a.size(), nb = (na + 1) >> 1;
    vector<int> b(a.begin(), a.begin() + nb);
    b = inv(b);
    vector<int> c = fft::multiply_mod(b, b, MOD);
    c.resize(na);
    c = fft::multiply_mod(a, c, MOD);
    b.resize(na), c.resize(na);
    for(int i = 0; i < na; i++) {
      c[i] = (((2ll * b[i] - c[i]) % MOD) + MOD) % MOD;
    }
    return c;
  }
  
  // A = B * C + D (mod x^n) (n = A.size())
  // always use with the next function mod
  // make sure A.size() >= B.size() or else it will return an empty vector
  vector<int> divide(const vector<int> &a, const vector<int> &b) {
    const int n = a.size(), m = b.size();
    if(n < m) return {};
    vector<int> A(a), B(b);
    reverse(A.begin(), A.end()), reverse(B.begin(), B.end());
    A.resize(n - m + 1), B.resize(n - m + 1);
    B = inv(B);
    vector<int> C = fft::multiply_mod(A, B, MOD);
    C.resize(n - m + 1), reverse(C.begin(), C.end());
    return C;
  }
  
  vector<int> mod(const vector<int> &a, const vector<int> &b, const vector<int> &c) {
    const int n = a.size(), m = b.size();
    if(n < m) return a;
    vector<int> e = fft::multiply_mod(b, c, MOD);
    e.resize(m - 1);
    for(int i = 0; i < m - 1; i++) {
      e[i] = a[i] - e[i];
      if(e[i] < 0) {
        e[i] += MOD;
      }
    }
    return e;
  }

  vector<int> buildPoly(const vector<int> &vec, const int left, const int right) {
    if(left == right) {
      vector<int> ret;
      ret.push_back(MOD - vec[left]);
      ret.push_back(1);
      return ret;
    }
    const int mid = (left + right) >> 1;
    return fft::multiply_mod(buildPoly(vec, left, mid), buildPoly(vec, mid + 1, right), MOD);
  }

  void multipointCalc(const vector<int> &poly, const vector<int> &vec, const int left, const int right, vector<int> &ret) {
    const int n = poly.size(), mid = (left + right) >> 1;
    if(n == 1) {
      for(int i = left; i <= right; i++) {
        ret[i] = poly[0];
      }
      return;
    }
    const vector<int> b0 = buildPoly(vec, left, mid);
    multipointCalc(mod(poly, b0, divide(poly, b0)), vec, left, mid, ret);
    if(left != right) {
      const vector<int> b1 = buildPoly(vec, mid + 1, right);
      multipointCalc(mod(poly, b1, divide(poly, b1)), vec, mid + 1, right, ret);
    }
  }

  // 多点求值
  vector<int> multipointCalc(const vector<int> &poly, const vector<int> &vec) {
    const int n = vec.size();
    vector<int> ret(n);
    multipointCalc(poly, vec, 0, n - 1, ret);
    return ret;
  }

  vector<int> multiInv(const vector<int> &vec) {
    const int n = vec.size();
    vector<int> a(n + 1), ret(n); a[0] = 1;
    for(int i = 1; i <= n; i++) {
      a[i] = 1ll * a[i - 1] * vec[i - 1] % MOD;
    }
    int cur = (exgcd(a[n], MOD).first + MOD) % MOD;
    for(int i = n - 1; i >= 0; i--) {
      ret[i] = 1ll * cur * a[i] % MOD;
      cur = 1ll * cur * vec[i] % MOD;
    }
    return ret;
  }

  // 快速插值 {{x0, y0}, {x1, y1}, {x2, y2}, ...}
  vector<int> interpolate(const vector<pair<int, int>> &p) {
    const int n = p.size(), n0 = (n + 1) >> 1, n1 = n - n0;
    if(n == 1) {
      return {p[0].second};
    }
    vector<pair<int, int>> p0(p.begin() + n1, p.end());
    vector<int> f0 = interpolate(p0);
    vector<int> x(n);
    for(int i = 0; i < n; i++) {
      x[i] = p[i].first;
    }
    vector<int> g0 = buildPoly(x, n1, n - 1);
    x.resize(n1);
    vector<int> fx = multipointCalc(f0, x), gx = multipointCalc(g0, x);
    gx = multiInv(gx);
    p0.resize(n1);
    for(int i = 0; i < n1; i++) {
      p0[i].first = p[0].first;
      p0[i].second = (p[i].second - fx[i] + MOD) % MOD;
      p0[i].second = 1ll * p0[i].second * gx[i] % MOD;
    }
    fx = interpolate(p0);
    fx = fft::multiply_mod(fx, g0, MOD);
    fx.resize(n), f0.resize(n);
    for(int i = 0; i < n; i++) {
      fx[i] = (fx[i] + f0[i]) % MOD;
    }
    return fx;
  }
}
