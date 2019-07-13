// reference exgcd and fft::multiply_mod
// make sure your nums are values in [0, MOD)
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
  // make sure A.size() >= B.size() or else if will return an empty vector
  vector<int> divide(const vector<int> &a, const vector<int> &b) {
    const int n = a.size(), m = b.size();
    if(n < m) {
      return {};
    }
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

  // 求导
  vector<int> derivate(const vector<int> &vec) {
    const int n = vec.size();
    if(n <= 1) {
      return {0};
    }
    vector<int> v(n - 1);
    for(int i = 1; i < n; i++) {
      v[i - 1] = 1ll * vec[i] * i % MOD;
    }
    return v;
  }

  // 积分
  // vector<int> intergral(const vector<int> &vec) {

  // }

  // 开根号
  // vector<int> sqrt(const vector<int> &vec) {
    
  // }

  // 对数函数 指数函数 快速幂

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

  // 快速插值 <x, y>
  // vector<int> interpolate(const vector<pair<int, int>> &p) {
  //   const int n = p.size();
  //   vector<int> x(n);
  //   for(int i = 0; i < n; i++) {
  //     x[i] = p[i].first;
  //   }
  //   vector<int> fm = buildPoly(x, 0, n - 1);
  //   fm = derivate(fm);
  //   fm = multipointCalc(fm, x);
  //   // above is for 分母
  //   // ...
  // }
}
