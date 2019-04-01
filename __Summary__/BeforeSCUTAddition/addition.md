## 线性筛（素数表、欧拉函数、莫比乌斯函数）

```cpp
    struct Seive {
        int maxn;
        vector<bool> isp;
        vector<int> p, phi, mu;

        Seive(int n = 0) : maxn(n), isp(n + 5, true), phi(n + 5, 0), mu(n + 5, 0) { solve(); }

        void solve() {
            isp[0] = isp[1] = false; phi[1] = 1; mu[1] = 1;
            for (int i = 2; i <= maxn; i++) {
                if (isp[i]) {
                    p.push_back(i);
                    phi[i] = i - 1;
                    mu[i] = -1;
                }
                for (int j = 0; j < (int)p.size() && i * p[j] <= maxn; j++) {
                    const int cur = i * p[j];
                    isp[cur] = false;
                    if (i % p[j]) {
                        phi[cur] = phi[i] * (p[j] - 1);
                        mu[cur] = -mu[i];
                    }
                    else {
                        phi[cur] = phi[i] * p[j];
                        mu[cur] = 0;
                        break;
                    }
                }
            }
        }
    };
```

## 扩展Lucas求大组合数

```cpp
    ll C(ll n, ll m, ll p) {
        if(m > n) return 0;
        ll ret = 1;
        for(ll i = 1; i <= m; i++) {
            ll a = (n + 1 - i) % p, b = mod(exgcd(i % p, p).first, p);
            ret = ret * a % p * b % p;
        }
        return ret;
    }

    ll lucas(ll n, ll m, ll p) {
        if(m == 0) {
            return 1;
        }
        return lucas(n / p, m / p, p) * C(n % p, m % p, p) % p;
    }

    ll cal(ll n, ll a, ll b, ll p) {
        if(!n) return 1;
        ll y = n / p, tmp = 1;
        for(ll i = 1; i <= p; i++) {
            if(i % a) {
                tmp = tmp * i % p;
            }
        }
        ll ans = fpow(tmp, y, p);
        for(ll i = y * p + 1; i <= n; i++) {
            if(i % a) {
                ans = ans * (i % p) % p;
            }
        }
        return ans * cal(n / a, a, b, p) % p;
    }

    ll multilucas(ll n, ll m, ll a, ll b, ll p) {
        ll s = 0;
        for(ll i = n; i; i /= a) s += i / a;
        for(ll i = m; i; i /= a) s -= i / a;
        for(ll i = n - m; i; i /= a) s -= i / a;
        ll tmp = fpow(a, s, p);
        ll t1 = cal(n, a, b, p), t2 = cal(m, a, b, p), t3 = cal(n - m, a, b, p);
        return tmp * t1 % p * mod(exgcd(t2, p).first, p) % p * mod(exgcd(t3, p).first, p) % p;
    }

    ll exlucas(ll n, ll m, ll p) {
        vector<ll>q, a;
        for(ll i = 2; i * i <= p; i++) {
            if(p % i == 0) {
                q.push_back(1);
                ll t = 0;
                while(p % i == 0) {
                    p /= i;
                    q.back() *= i;
                    t++;
                }
                a.push_back(q.back() == i ? lucas(n, m, q.back()) : multilucas(n, m, i, t, q.back()));
            }
        }
        if(p > 1) {
            q.push_back(p);
            a.push_back(lucas(n, m, p));
        }
        const int e = q.size();
        for(ll i = 1; i < e; i++) {
            pll d = exgcd(q[0], q[i]);
            ll c = a[i] - a[0], g = d.first * q[0] + d.second * q[i];
            if(c % g) exit(-1);
            a[0] = q[0] * mod(c / g * d.first, q[i] / g) + a[0];
            q[0] = q[0] * q[i] / g;
        }
        return a[0];
    }
```

## exBSGS

```cpp
    ll bsgs(ll a, ll b, ll c, ll q = 1, ll d = 0) {
        unordered_map<ll, ll> x;
        ll m = sqrt(c) + 1;
        ll v = 1;
        if(d > 0) {
            for(int i = 1; i <= m; i++) {
                v = fmul(v, a, c);
                x[fmul(v, b, c)] = i;
            }
        } else {
            for(int i = 0; i < m; i++) {
                x[fmul(v, b, c)] = i;
                v = fmul(v, a, c);
            }
        }
        for(int i = 1; i <= m; i++) {
            q = fmul(q, v, c);
            auto it = x.find(q);
            if(it != x.end()) {
                return i * m - it->second + d;
            }
        }
        return -1;
    }

    ll exbsgs(ll a, ll b, ll m) {
        a = mod(a, m), b = mod(b, m);
        if(a == 0) {
            return b > 1 ? -1 : b == 0 && m > 1;
        }
        if(b == 1 && gcd(a, m) != 1) { // b为1时随机应变吧。
            return -1;
        }
        ll g, c = 0, q = 1;
        while((g = gcd(a, m)) != 1) {
            if(b == 1) return c;
            if(b % g) return -1;
            c++;
            b /= g, m /= g;
            q = fmul(a / g, q, m);
        }
        return bsgs(a, b, m, q, c);
    }
```

## 任意模数fft by tourist

```cpp
namespace fft {
  const double pi = acos(-1.0);
  struct Complex {
    double r, i;
    Complex(double x = 0, double y = 0) : r(x), i(y) {}
    Complex operator+ (const Complex& b) const {
      return Complex(r + b.r, i + b.i);
    }
    Complex operator- (const Complex& b) const {
      return Complex(r - b.r, i - b.i);
    }
    Complex operator* (const Complex& b) const {
      return Complex(r * b.r - i * b.i, r * b.i + i * b.r);
    }
  };
  Complex conj(Complex a) { return Complex(a.r, -a.i); }

  int base = 1;
  vector<int> rev = { 0, 1 };
  vector<Complex> roots = { { 0, 0 }, { 1, 0 } };
  
  void ensure_base(int nbase) {
    if (nbase <= base) return;
    rev.resize(1 << nbase);
    for (int i = 0; i < (1 << nbase); i++) {
      rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
    }
    roots.resize(1 << nbase);
    while (base < nbase) {
      double angle = 2 * pi / (1 << (base + 1));
      for (int i = 1 << (base - 1); i < (1 << base); i++) {
        roots[i << 1] = roots[i];
        double angle_i = angle * (2 * i + 1 - (1 << base));
        roots[(i << 1) + 1] = Complex(cos(angle_i), sin(angle_i));
      }
      base++;
    }
  }
  void fft(vector<Complex> &a, int n = -1) {
    if (n == -1) {
      n = a.size();
    }
    assert((n & (n - 1)) == 0);
    int zeros = __builtin_ctz(n);
    ensure_base(zeros);
    int shift = base - zeros;
    for (int i = 0; i < n; i++) {
      if (i < (rev[i] >> shift)) {
        swap(a[i], a[rev[i] >> shift]);
      }
    }
    for (int k = 1; k < n; k <<= 1) {
      for (int i = 0; i < n; i += 2 * k) {
        for (int j = 0; j < k; j++) {
          Complex z = a[i + j + k] * roots[j + k];
          a[i + j + k] = a[i + j] - z;
          a[i + j] = a[i + j] + z;
        }
      }
    }
  }
  vector<Complex> fa, fb;
  vector<int> multiply(vector<int> &a, vector<int> &b) {
    int need = a.size() + b.size() - 1;
    int nbase = 32 - __builtin_clz(need);
    ensure_base(nbase);
    int sz = 1 << nbase;
    if (sz > (int) fa.size()) {
      fa.resize(sz);
    }
    for (int i = 0; i < sz; i++) {
      int x = (i < (int) a.size() ? a[i] : 0);
      int y = (i < (int) b.size() ? b[i] : 0);
      fa[i] = Complex(x, y);
    }
    fft(fa, sz);
    Complex r(0, -0.25 / sz);
    for (int i = 0; i <= (sz >> 1); i++) {
      int j = (sz - i) & (sz - 1);
      Complex z = (fa[j] * fa[j] - conj(fa[i] * fa[i])) * r;
      if (i != j) {
        fa[j] = (fa[i] * fa[i] - conj(fa[j] * fa[j])) * r;
      }
      fa[i] = z;
    }
    fft(fa, sz);
    vector<int> res(need);
    for (int i = 0; i < need; i++) {
      res[i] = fa[i].r + 0.5;
    }
    return res;
  }
  vector<int> multiply_mod(vector<int> &a, vector<int> &b, int m, int eq = 0) {
    int need = a.size() + b.size() - 1;
    int nbase = 31 - __builtin_clz(need);
    ensure_base(nbase);
    int sz = 1 << nbase;
    if (sz > (int) fa.size()) {
      fa.resize(sz);
    }
    for (int i = 0; i < (int) a.size(); i++) {
      int x = (a[i] % m + m) % m;
      fa[i] = Complex(x & ((1 << 15) - 1), x >> 15);
    }
    fill(fa.begin() + a.size(), fa.begin() + sz, Complex {0, 0});
    fft(fa, sz);
    if (sz > (int) fb.size()) {
      fb.resize(sz);
    }
    if (eq) {
      copy(fa.begin(), fa.begin() + sz, fb.begin());
    } else {
      for (int i = 0; i < (int) b.size(); i++) {
        int x = (b[i] % m + m) % m;
        fb[i] = Complex(x & ((1 << 15) - 1), x >> 15);
      }
      fill(fb.begin() + b.size(), fb.begin() + sz, Complex {0, 0});
      fft(fb, sz);
    }
    double ratio = 0.25 / sz;
    Complex r2(0, -1), r3(ratio, 0), r4(0, -ratio), r5(0, 1);
    for (int i = 0; i <= (sz >> 1); i++) {
      int j = (sz - i) & (sz - 1);
      Complex a1 = (fa[i] + conj(fa[j]));
      Complex a2 = (fa[i] - conj(fa[j])) * r2;
      Complex b1 = (fb[i] + conj(fb[j])) * r3;
      Complex b2 = (fb[i] - conj(fb[j])) * r4;
      if (i != j) {
        Complex c1 = (fa[j] + conj(fa[i]));
        Complex c2 = (fa[j] - conj(fa[i])) * r2;
        Complex d1 = (fb[j] + conj(fb[i])) * r3;
        Complex d2 = (fb[j] - conj(fb[i])) * r4;
        fa[i] = c1 * d1 + c2 * d2 * r5;
        fb[i] = c1 * d2 + c2 * d1;
      }
      fa[j] = a1 * b1 + a2 * b2 * r5;
      fb[j] = a1 * b2 + a2 * b1;
    }
    fft(fa, sz);
    fft(fb, sz);
    vector<int> res(need);
    for (int i = 0; i < need; i++) {
      long long aa = fa[i].r + 0.5;
      long long bb = fb[i].r + 0.5;
      long long cc = fa[i].i + 0.5;
      res[i] = (aa + ((bb % m) << 15) + ((cc % m) << 30)) % m;
    }
    return res;
  }
  vector<int> square_mod(vector<int> &a, int m) {
    return multiply_mod(a, a, m, 1);
  }
}
```

## ntt

```cpp
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
    int sz = (int)a.size() + (int)b.size() + 1;
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
```

## 区间增区间查bit

```cpp
// ~i & (i + 1)
struct Interval {
  int N, base[2][maxn];
  void setN(int n) { N = n; }
  void init() { memset(base, 0, sizeof base); }
  void add(int at, int v) {
    if (!at) return;
    for (int i = at; i <= N; i += i & -i) {
      base[0][i] += v, base[1][i] -= v * at;
    }
  }
  void add(int l, int r, int v) {
    add(l, v), add(r + 1, -v);
  }
  int getSum(int at) {
    int sum = 0, mul = at + 1;
    for (int i = at; i; i -= i & -i) {
      sum += mul * base[0][i] + base[1][i];
    }
    return sum;
  }
  int query(int l, int r) {
    return getSum(r) - getSum(l - 1);
  }
};
```

## Suffix Automa

```cpp
struct SAM {
  int last, tot, sz[maxn << 1], len[maxn << 1], fa[maxn << 1];
  int ch[maxn << 1][30];

  SAM() {
    tot = 0, last = newNode(0), len[0] = -1;
    memset(sz, 0, sizeof sz);
  }
  inline int newNode(int v) {
    len[++tot] = v, fa[tot] = 0;
    memset(ch[tot], 0, sizeof ch[tot]);
    return tot;
  }
  void append(int c) {
    int p = last, u = newNode(len[last] + 1);
    for (; p && !ch[p][c]; p = fa[p]) {
      ch[p][c] = u;
    }
    if (p == 0) {
      fa[u] = 1;
    } else {
      int q = ch[p][c];
      if (len[q] == len[p] + 1) {
        fa[u] = q;
      } else {
        int nq = newNode(len[p] + 1);
        memcpy(ch[nq], ch[q], sizeof ch[q]);
        fa[nq] = fa[q], fa[u] = fa[q] = nq;
        for (; p && (ch[p][c] == q); p = fa[p]) {
          ch[p][c] = nq;
        }
      }
    }
    last = u;
  }
  void match(char *s) {
    int pos = 1, length = 0;
    for (int i = 0, n = strlen(s); i < n; i++) {
      while (pos && !ch[pos][s[i] - 'a']) {
        pos = fa[pos], length = len[pos];
      }
      if (pos) {
        ++length, pos = ch[pos][s[i] - 'a'];
        // update ans
      } else {
        pos = 1, length = 0;
      }
    }
  }
} sam;
```

## Suffix Array by tourist

```cpp
template <typename T>
vector<int> build_sa(int n, const T &s, int charset) {
  vector<int> a(n);
  if (n == 0) {
    return a;
  }
  if (charset != -1) {
    vector<int> aux(charset, 0);
    for (int i = 0; i < n; i++) {
      aux[s[i]]++;
    }
    int sum = 0;
    for (int i = 0; i < charset; i++) {
      int add = aux[i];
      aux[i] = sum;
      sum += add;
    }
    for (int i = 0; i < n; i++) {
      a[aux[s[i]]++] = i;
    }
  } else {
    iota(a.begin(), a.end(), 0);
    sort(a.begin(), a.end(), [&s](int i, int j) { return s[i] < s[j]; });
  }
  vector<int> sorted_by_second(n), ptr_group(n);
  vector<int> new_group(n), group(n);
  group[a[0]] = 0;
  for (int i = 1; i < n; i++) {
    group[a[i]] = group[a[i - 1]] + (!(s[a[i]] == s[a[i - 1]]));
  }
  int cnt = group[a[n - 1]] + 1;
  int step = 1;
  while (cnt < n) {
    int at = 0;
    for (int i = n - step; i < n; i++) {
      sorted_by_second[at++] = i;
    }
    for (int i = 0; i < n; i++) {
      if (a[i] - step >= 0) {
        sorted_by_second[at++] = a[i] - step;
      }
    }
    for (int i = n - 1; i >= 0; i--) {
      ptr_group[group[a[i]]] = i;
    }
    for (int i = 0; i < n; i++) {
      int x = sorted_by_second[i];
      a[ptr_group[group[x]]++] = x;
    }
    new_group[a[0]] = 0;
    for (int i = 1; i < n; i++) {
      if (group[a[i]] != group[a[i - 1]]) {
        new_group[a[i]] = new_group[a[i - 1]] + 1;
      } else {
        int pre = ((a[i - 1] + step >= n) ? -1 : group[a[i - 1] + step]);
        int cur = ((a[i] + step >= n) ? -1 : group[a[i] + step]);
        new_group[a[i]] = new_group[a[i - 1]] + (pre != cur);
      }
    }
    swap(group, new_group);
    cnt = group[a[n - 1]] + 1;
    step <<= 1;
  }
  return a;
}

template <typename T>
vector<int> build_sa(const T &s, int charset) {
  return build_sa((int)s.size(), s, charset);
}

template <typename T>
vector<int> build_lcp(int n, const T &s, const vector<int> &sa) {
  assert((int) sa.size() == n);
  vector<int> pos(n);
  for (int i = 0; i < n; i++) {
    pos[sa[i]] = i;
  }
  vector<int> lcp(max(n - 1, 0));
  int k = 0;
  for (int i = 0; i < n; i++) {
    k = max(k - 1, 0);
    if (pos[i] == n - 1) {
      k = 0;
    } else {
      int j = sa[pos[i] + 1];
      while (i + k < n && j + k < n && s[i + k] == s[j + k]) {
        k++;
      }
      lcp[pos[i]] = k;
    }
  }
  return lcp;
}

template <typename T>
vector<int> build_lcp(const T &s, const vector<int> &sa) {
  return build_lcp((int)s.size(), s, sa);
}
```

## Z-func

```cpp
// 与exkmp类似
// z[i]: suf_s[i..n-1] 与 s[0...n-1] 的公共前缀长度
vector<int> getZ(string s) {
  int n = (int)s.length();
  vector<int> z(n);
  for (int i = 1, l = 0, r = 0; i < n; i++) {
    if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) ++z[i];
    if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
  }
  return z;
}
```

## 2-sat

```cpp
struct twoSat {
  struct edge {
    int v, next;
    edge(int a = 0, int b = 0) : v(a), next(b) {}
  }G[maxm];
  int tot, head[maxn], mark[maxn], sz, stk[maxn];
  void init() {
    tot = 0;
    memset(mark, 0, sizeof mark);
    memset(head, -1, sizeof head);
  }
  // for every case u, (status[u] xor status[u ^ 1]) == true.
  // 
  // addcase: if status[u] == true then status[v] == true,
  // but if status[u] == false then status[v] can be true or false.
  // 
  void addcase(int u, int v) {
    G[tot] = edge(v, head[u]); head[u] = tot++;
  }
  int dfs(int u) {
    if (mark[u ^ 1]) return 0;
    if (mark[u]) return 1;
    stk[sz++] = u, mark[u] = 1;
    for (int i = head[u]; ~i; i = G[i].next) {
      if (!dfs(G[i].v)) return 0;
    }
    return 1;
  }
  int solve(int n) {
    for (int i = 0; i < n; i += 2) {
      if (!mark[i] && !mark[i ^ 1]) {
        sz = 0;
        if (!dfs(i)) {
          while (sz > 0) mark[stk[--sz]] = 0;
          if (!dfs(i ^ 1)) return 0;
        }
      }
    }
    return 1;
  }
}sat;
```

## 点双

```cpp
// vertex-bcc
struct bcc {
  struct edge { int u, v; };
  vector<int> G[maxn], cont[maxn];
  int N, tag, tot, dfn[maxn], bccno[maxn];
  bool iscut[maxn];
  stack<edge> S;

  void init(int n) {
    N = n, tag = tot = 0;
    for (int i = 1; i <= N; i++) {
      G[i].clear();
      dfn[i] = bccno[i] = 0;
      iscut[i] = false;
    }
  }
  void addedge(int u, int v) {
    G[u].push_back(v), G[v].push_back(u);
  }
  int dfs(int u, int f) {
    int lowu = dfn[u] = ++tag;
    int child = 0;
    for (auto& v : G[u]) {
      if (!dfn[v]) {
        ++child, S.push({ u, v });
        int lowv = dfs(v, u);
        lowu = min(lowu, lowv);
        if (lowv >= dfn[u]) {
          iscut[u] = true;
          cont[++tot].clear();
          while (true) {
            edge e = S.top(); S.pop();
            if (bccno[e.u] != tot) {
              cont[tot].push_back(e.u);
              bccno[e.u] = tot;
            }
            if (bccno[e.v] != tot) {
              cont[tot].push_back(e.v);
              bccno[e.v] = tot;
            }
            if (e.u == u && e.v == v) {
              break;
            }
          }
        }
      } else if (dfn[v] < dfn[u] && v != f) {
        S.push({ u, v });
        lowu = min(lowu, dfn[v]);
      }
    }
    if (f < 0 && child == 1) {
      iscut[u] = false;
    }
    return lowu;
  }
};
```

## 桥

```cpp
struct biBridge {
  vector<pair<int, int>> G[maxn];
  int N, tag, dfn[maxn], low[maxn];
  bool isbridge[maxm];

  void init(int n) {
    N = n, tag = 0;
    for (int i = 1; i <= n; i++) {
      dfn[i] = low[i] = 0, G[i].clear();
    }
    memset(isbridge, 0, sizeof isbridge);
  }
  void addedge(int u, int v, int i) {
    G[u].push_back({ v, i });
    G[v].push_back({ u, i });
  }
  void dfs(int u, int f) {
    dfn[u] = low[u] = ++tag;
    int child = 0;
    for (auto& e : G[u]) {
      int v = e.first, index = e.second;
      if (!dfn[v]) {
        ++child, dfs(v, u);
        low[u] = min(low[u], low[v]);
        if (low[v] > dfn[u]) {
          isbridge[index] = true;
        }
      } else if (dfn[v] < dfn[u] && v != f) {
        low[u] = min(low[u], dfn[v]);
      }
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!dfn[i]) dfs(i, -1);
    }
  }
} bridge;
```

## 割

```cpp
struct biCut {
  vector<int> G[maxn];
  int N, tag, dfn[maxn], low[maxn];
  bool iscut[maxn];

  void init(int n) {
    N = n, tag = 0;
    for (int i = 1; i <= N; i++) {
      dfn[i] = low[i] = 0, iscut[i] = false;
      G[i].clear();
    }
  }
  void addedge(int u, int v) {
    G[u].push_back(v), G[v].push_back(u);
  }
  void dfs(int u, int f) {
    low[u] = dfn[u] = ++tag;
    int child = 0;
    for (auto& v : G[u]) {
      if (!dfn[v]) {
        ++child, dfs(v, u);
        low[u] = min(low[u], low[v]);
        if (low[v] >= dfn[u]) {
          iscut[u] = true;
        }
      } else if (dfn[v] < dfn[u] && v != f) {
        low[u] = min(low[u], dfn[v]);
      }
    }
    if (f < 0 && child == 1) {
      iscut[u] = false;
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!dfn[i]) dfs(i, -1);
    }
  }
};
```

## Kosaraju

```cpp
struct kosaraju {
  int N, tot, scc[maxn], vis[maxn];
  vector<int> G[maxn], R[maxn], acc;
  void init(int n) {
    N = n;
    tot = 0, acc.clear();
    for (int i = 1; i <= N; i++) {
      G[i].clear(), R[i].clear();
      vis[i] = 0, scc[i] = 0;
    }
  }
  void DFS1(int u) {
    vis[u] = 1;
    for (auto& v : G[u]) {
      if (!vis[v]) DFS1(v);
    }
    acc.push_back(u);
  }
  void DFS2(int u, int p) {
    scc[u] = p;
    for (auto& v : R[u]) {
      if (!scc[v]) DFS2(v, p);
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!vis[i]) DFS1(i);
    }
    reverse(acc.begin(), acc.end());
    for (auto& u : acc) {
      if (!scc[u]) DFS2(u, ++tot);
    }
  }
};
```

