typedef long long ll;
typedef vector<ll> VL;
typedef pair<ll, ll> pll;
namespace crt {
  pll exgcd(const ll x, const ll y) {
    if(!y) return {1, 0};
    pll cur = exgcd(y, x % y);
    return {cur.second, cur.first - (x / y) * cur.second};
  }
  ll mod(ll x, ll p) {
    x %= p;
    return x >= 0 ? x : x + p;
  }
  pll crt(const vector<pll> &v) {
    ll a = 1, r = 0;
    const int n = v.size();
    for(int i = 0; i < n; i++) {
      pll cur = exgcd(a, v[i].first);
      ll g = a * cur.first + v[i].first * cur.second;
      if((v[i].second - r) % g != 0) return {-1, -1};
      const ll p = v[i].first / g;
      r += mod(cur.first * ((v[i].second - r) / g), p) * a;
      a *= p;
    }
    return {a, r};
  }
}

namespace LinearRecurrence {
  void extand(VL &v, size_t d, ll value = 0) {
    if(d <= v.size()) return;
    v.resize(d, value);
  }
  VL ReedsSloane(const VL &s, ll Mod) {
    function<ll(ll, ll)> inverse = [](ll a, ll m) {
      pll ret = crt::exgcd(a, m);
      ll g = a * ret.first + m * ret.second;
      if(g != 1) return -1ll;
      return ret.first >= 0 ? ret.first : ret.first + m;
    };
    function<int(const VL&, const VL&)> L = [](const VL &a, const VL &b) {
      int da = (a.size() > 1 || (a.size() == 1 && a[0])) ? a.size() - 1 : -1000;
      int db = (b.size() > 1 || (b.size() == 1 && b[0])) ? b.size() - 1 : -1000;
      return max(da, db + 1);
    };
    function<pair<VL, VL>(const VL&, ll, ll, ll)> prime_power = [&](const VL &s, ll Mod, ll p, ll e) {
      vector<VL> a(e), b(e), an(e), bn(e), ao(e), bo(e);
      VL t(e), u(e), r(e), to(e, 1), uo(e), pw(e + 1);
      pw[0] = 1;
      for(int i = pw[0] = 1; i <= e; i++) pw[i] = pw[i - 1] * p;
      for(ll i = 0; i < e; i++) {
        a[i] = {pw[i]}, an[i] = {pw[i]};
        b[i] = {0}, bn[i] = {s[0] * pw[i] % Mod};
        t[i] = s[0] * pw[i] % Mod;
        if(t[i] == 0) t[i] = 1, u[i] = e;
        else {
          for(u[i] = 0; t[i] % p == 0; t[i] /= p, u[i]++);
        }
      }
      for(size_t k = 1; k < s.size(); k++) {
        for(int g = 0; g < e; g++) {
          if(L(an[g], bn[g]) > L(a[g], b[g])) {
            ao[g] = a[e - 1 - u[g]];
            bo[g] = b[e - 1 - u[g]];
            to[g] = t[e - 1 - u[g]];
            uo[g] = u[e - 1 - u[g]];
            r[g] = k - 1;
          }
        }
        a = an, b = bn;
        for(int o = 0; o < e; o++) {
          ll d = 0;
          for(size_t i = 0; i < a[o].size() && i <= k; i++) {
            d = (d + a[o][i] * s[k - i]) % Mod;
          }
          if(d == 0) t[o] = 1, u[o] = e;
          else {
            for(u[o] = 0, t[o] = d; t[o] % p == 0; t[o] /= p, u[o]++);
            int g = e - 1 - u[o];
            if(L(a[g], b[g]) == 0) {
              extand(bn[o], k + 1);
              bn[o][k] = (bn[o][k] + d) % Mod;
            } else {
              ll coef = t[o] * inverse(to[g], Mod) % Mod * pw[u[o] - uo[g]] % Mod;
              int m = k - r[g];
              extand(an[o], ao[g].size() + m);
              extand(bn[o], bo[g].size() + m);
              for(size_t i = 0; i < ao[g].size(); i++) {
                an[o][i + m] -= coef * ao[g][i] % Mod;
                if(an[o][i + m] < 0) an[o][i + m] += Mod;
              }
              while(an[o].size() && an[o].back() == 0) an[o].pop_back();
              for(size_t i = 0; i < bo[g].size(); i++) {
                bn[o][i + m] -= coef * bo[g][i] % Mod;
                if(bo[o][i + m] < 0) bn[o][i + m] += Mod;
              }
              while(bn[o].size() && bn[o].back() == 0) bn[o].pop_back();
            }
          }
        }
      }
      return make_pair(an[0], bn[0]);
    };

    vector<tuple<ll, ll, int>> fac;
    for(ll i = 2; i * i <= Mod; i++) {
      if(Mod % i == 0) {
        ll cnt = 0, pw = 1;
        while(Mod % i == 0) Mod /= i, ++cnt, pw *= i;
        fac.emplace_back(pw, i, cnt);
      }
    }
    if(Mod > 1) fac.emplace_back(Mod, Mod, 1);
    vector<VL> as;
    size_t n = 0;
    for(auto &&x: fac) {
      ll mod, p, e;
      VL a, b;
      tie(mod, p, e) = x;
      auto ss = s;
      for(auto &&x: ss) x %= mod;
      tie(a, b) = prime_power(ss, mod, p, e);
      as.emplace_back(a);
      n = max(n, a.size());
    }
    VL a(n);
    vector<pll> c(as.size());
    for(size_t i = 0; i < n; i++) {
      for(size_t j = 0; j < as.size(); j++) {
        c[j].first = get<0>(fac[j]);
        c[j].second = i < as[j].size() ? as[j][i] : 0;
      }
      a[i] = crt::crt(c).second;
    }
    return a;
  }

  int m;
  VL ini, trans;
  ll Mod;
  void init(const VL &s, ll mod) {
    Mod = mod;
    VL A = ReedsSloane(s, Mod);
    if(A.empty()) A = {0};
    m = A.size() - 1;
    trans.resize(m);
    for(int i = 0; i < m; i++) {
      trans[i] = (Mod - A[i + 1]) % Mod;
    }
    reverse(trans.begin(), trans.end());
    ini = {s.begin(), s.begin() + m};
  }

  ll calc(ll n) {
    if(Mod == 1) return n;
    if(n < m) return ini[n];
    VL v(m), u(m << 1);
    int msk = !!n;
    for(ll m = n; m > 1; m >>= 1) msk <<= 1;
    v[0] = 1 % Mod;
    for(int x = 0; msk; msk >>= 1, x <<= 1) {
      fill_n(u.begin(), m * 2, 0);
      x |= !!(n & msk);
      if(x < m) u[x] = 1 % Mod;
      else {
        for(int i = 0; i < m; i++) {
          for(int j = 0, t = i + (x & 1); j < m; j++, t++) {
            u[t] = (u[t] + v[i] * v[j]) % Mod; // to better
          }
        }
        for(int i = m * 2 - 1; i >= m; i--) {
          for(int j = 0, t = i - m; j < m; j++, t++) {
            u[t] = (u[t] + trans[j] * u[i]) % Mod; // to better
          }
        }
      }
      v = {u.begin(), u.begin() + m};
    }
    ll ret = 0;
    for(int i = 0; i < m; i++) {
      ret = (ret + v[i] * ini[i]) % Mod; // to better
    }
    return ret;
  }
}