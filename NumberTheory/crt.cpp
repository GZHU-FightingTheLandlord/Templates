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
    //v里每个pll中first为被模数，second为模数
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