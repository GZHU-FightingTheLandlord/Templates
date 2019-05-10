## 数学与数论

### 自适应Simpson积分

$\int_a^bF(x)dx\Rightarrow$ `asr(a, b, eps, simpson(a, b))`

```cpp
double simpson(const double& a, const double& b) {
  double c = (a + b) / 2;
  return (F(a) + 4 * F(c) + F(b)) * (b - a) / 6;
}
double asr(double a, double b, double eps, double A) {
  double c = (a + b) / 2;
  double L = simpson(a, c), R = simpson(c, b);
  if (fabs(L + R - A) <= 15 * eps)
    return L + R + (L + R - A) / 15.0;
  return asr(a, c, eps / 2, L) + asr(c, b, eps / 2, R);
}
```

### BM推公式大法

```cpp
struct BM {
  static const int MAXN = 10005;
  int n, pn, fail[MAXN];
  double delta[MAXN];
  vector<double> ps[MAXN];
  void Solve(double x[], const int &n) {
    pn = 0;
    memset(fail, 0, sizeof fail);
    memset(delta, 0, sizeof delta);
    ps[0].clear();
    for (int i = 1; i <= n; i++) {
      double dt = -x[i];
      for (int j = 0; j < ps[pn].size(); j++) {
        dt += x[i - j - 1] * ps[pn][j];
      }
      delta[i] = dt;
      if (fabs(dt) <= 1e-8) continue;
      fail[pn] = i;
      if (!pn) {
        ps[++pn].resize(1);
        continue;
      }
      vector<double> &ls = ps[pn - 1];
      double k = -dt / delta[fail[pn - 1]];
      vector<double> cur(i - fail[pn - 1] - 1);
      cur.push_back(-k);
      for (int j = 0; j < ls.size(); j++) {
        cur.push_back(ls[j] * k);
      }
      if (cur.size() < ps[pn].size()) {
        cur.resize(ps[pn].size());
      }
      for (int j = 0; j < ps[pn].size(); j++) {
        cur[j] += ps[pn][j];
      }
      ps[++pn] = cur;
    }
  }
  void print() {
    for (int i = 0; i < ps[pn].size(); i++) {
      printf("%lf%c", ps[pn][i], (i == ps[pn].size() - 1) ? '\n' : ' ');
    }
  }
} B;
double x[BM::MAXN];
int main() {
  for (int n; ~scanf("%d", &n); ) {
    for (int i = 1; i <= n; i++) {
      scanf("%lf", &x[i]);
    }
    B.Solve(x, n), B.print();
  }
}
```

### 线性基

```cpp
struct LinearBasis {
  const static int MAXL = 50;
  long long a[MAXL + 1];
  LinearBasis() {
    memset(a, 0, sizeof a);
  }
  void insert(long long t) {
    for (int j = MAXL; j >= 0; j--) {
      if (!(t & (1ll << j))) continue;
      if (a[j]) t ^= a[j];
      else {
        for (int k = 0; k < j; k++) if (t & (1ll << k)) t ^= a[k];
        for (int k = j + 1; k <= MAXL; k++) if (a[k] & (1ll << j)) a[k] ^= t;
        a[j] = t;
        return;
      }
    }
  }
};
```

### 扩展欧几里得

```cpp
pll exgcd(const long long x, const long long y) {
  if (!y) return make_pair(1, 0);
  pll cur = exgcd(y, x % y);
  return make_pair(cur.second, cur.first - (x / y) * cur.second);
}
```

### 中国剩余定理

```cpp
//v里每个pll中first为被模数，second为模数
pll crt(const vector<pll> & v) {
  ll a = 1, r = 0;
  const int len = v.size();
  for(int i = 0; i < len; i++) {
    pll cur = exgcd(a, v[i].first);
    ll gcd = a * cur.first + v[i].first * cur.second;
    if((v[i].second - r) % gcd != 0) {
      return make_pair(-1, -1);
    }
    const ll p = v[i].first / gcd;
    r += mod(cur.first * ((v[i].second - r) / gcd), p) * a;
    a *= p;
  }
  return make_pair(a, r);
}
```

### 扩展卢卡斯

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

### 快速乘

```cpp
// mod <= 1e12
inline ll mul(ll a, ll b, ll mod) {
  return (((a * (b >> 20) % mod) << 20) + (a * (b & ((1 << 20) - 1)))) % mod;
}
// mod <= 1e18
inline ll mul(ll a, ll b, ll mod) {
  ll d = (ll)floor(a * (long double)b / mod + 0.5);
  ll ret = (a * b - d * mod) % mod;
  if (ret < 0) ret += mod;
  return ret;
}
```

### exbsgs

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
// 返回最小正整数n使得 a^n mod m = b; O(sqrt(m))
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

### polysum

```cpp
namespace polysum {
ll mod = 998244353LL;
#define rep(i,a,n) for (int i=a;i<n;i++)
#define per(i,a,n) for (int i=n-1;i>=a;i--)
const int D=200005;
ll a[D],f[D],g[D],p[D],p1[D],p2[D],b[D],h[D][2],C[D];
ll powmod(ll a,ll b) {
  ll res=1;
  a%=mod;
  assert(b>=0);
  for(; b; b>>=1) {
    if(b&1)res=res*a%mod;
    a=a*a%mod;
  }
  return res;
}
//函数用途：给出数列的（d+1）项，其中d为最高次方项
//求出数列的第n项，数组下标从0开始
ll calcn(int d,ll *a,ll n) { // a[0].. a[d]  a[n]
  if (n<=d) return a[n];
  p1[0]=p2[0]=1;
  rep(i,0,d+1) {
    ll t=(n-i+mod)%mod;
    p1[i+1]=p1[i]*t%mod;
  }
  rep(i,0,d+1) {
    ll t=(n-d+i+mod)%mod;
    p2[i+1]=p2[i]*t%mod;
  }
  ll ans=0;
  rep(i,0,d+1) {
    ll t=g[i]*g[d-i]%mod*p1[i]%mod*p2[d-i]%mod*a[i]%mod;
    if ((d-i)&1) ans=(ans-t+mod)%mod;
    else ans=(ans+t)%mod;
  }
  return ans;
}
void init(int M) {
  f[0]=f[1]=g[0]=g[1]=1;
  rep(i,2,M+5) f[i]=f[i-1]*i%mod;
  g[M+4]=powmod(f[M+4],mod-2);
  per(i,1,M+4) g[i]=g[i+1]*(i+1)%mod;
}
//函数用途：给出数列的（m+1）项，其中m为最高次方
//求出数列的前（n-1）项的和
ll polysum(ll m,ll *a,ll n) { // a[0].. a[m] \sum_{i=0}^{n-1} a[i]
  ll b[D];
  for(int i=0; i<=m; i++) b[i]=a[i];
  b[m+1]=calcn(m,b,m+1);
  rep(i,1,m+2) b[i]=(b[i-1]+b[i])%mod;
  return calcn(m+1,b,n-1);
}
ll qpolysum(ll R,ll n,ll *a,ll m) { // a[0].. a[m] \sum_{i=0}^{n-1} a[i]*R^i
  if (R==1) return polysum(n,a,m);
  a[m+1]=calcn(m,a,m+1);
  ll r=powmod(R,mod-2),p3=0,p4=0,c,ans;
  h[0][0]=0;
  h[0][1]=1;
  rep(i,1,m+2) {
    h[i][0]=(h[i-1][0]+a[i-1])*r%mod;
    h[i][1]=h[i-1][1]*r%mod;
  }
  rep(i,0,m+2) {
    ll t=g[i]*g[m+1-i]%mod;
    if (i&1) p3=((p3-h[i][0]*t)%mod+mod)%mod,p4=((p4-h[i][1]*t)%mod+mod)%mod;
    else p3=(p3+h[i][0]*t)%mod,p4=(p4+h[i][1]*t)%mod;
  }
  c=powmod(p4,mod-2)*(mod-p3)%mod;
  rep(i,0,m+2) h[i][0]=(h[i][0]+h[i][1]*c)%mod;
  rep(i,0,m+2) C[i]=h[i][0];
  ans=(calcn(m,C,n)*powmod(R,n)-c)%mod;
  if (ans<0) ans+=mod;
  return ans;
}
} // polysum::init();
```


### 线性筛

```cpp
struct Seive {
  int maxn;
  vector<bool> isp;
  vector<int> p, phi, mu;
  Seive(int n = 0) : maxn(n), isp(n + 5, true), phi(n + 5, 0), mu(n + 5, 0) { solve(); }
  void solve() {
    isp[0] = isp[1] = false;
    phi[1] = 1;
    mu[1] = 1;
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
        } else {
          phi[cur] = phi[i] * p[j];
          mu[cur] = 0;
          break;
        }
      }
    }
  }
};
```

### MillerRabin素性测试

```cpp
const int psize = 1010000;
bool isp[psize];
int prime[psize], tot;
void prime_table() {
  register int i, j;
  for (i = 2, tot = 0; i < psize; i++) {
    if (!isp[i]) prime[tot++] = i;
    for (j = 0; j < tot && prime[j] * i < psize; j++) {
      isp[prime[j] * i] = true;
      if (i % prime[j] == 0) break;
    }
  }
}
bool witness(ll a, ll n) {
  int t = 0;
  ll u = n - 1;
  for (; ~u & 1; u >>= 1) t++;
  ll x = qpow(a, u, n), _x = 0;
  while (t--) {
    _x = mul(x, x, n);
    if (_x == 1 && x != 1 && x != n - 1) return true;
    x = _x;
  }
  return _x != 1;
}
bool Miller(ll n) {
  if (n < 2) return false;
  if (n < psize) return !isp[n];
  if (~n & 1) return false;
  for (int j = 0; j <= 7; j++) {
    if (witness(rand() % (n - 1) + 1, n)) {
      return false;
    }
  }
  return true;
}
```

### pollard_rho分解质因数

```cpp
int tot;
long long factor[10000];
long long pollard_rho(long long x, long long c) {
  long long i = 1, k = 2;
  long long x0 = rand() % x, y = x0;
  while (true) {
    i++;
    x0 = (mul(x0, x0, x) + c) % x;
    long long d = __gcd(y - x0, x);
    if (d != 1 && d != x) return d;
    if (y == x0) return x;
    if (i == k) {
      y = x0, k <<= 1;
    }
  }
}
void findfac(long long n) {
  if (Miller(n)) {
    factor[tot++] = n;
    return;
  }
  long long p = n;
  while (p >= n) p = pollard_rho(p, rand() % (n - 1) + 1);
  findfac(p), findfac(n / p);
}
```

### fft

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
  vector<int> multiply(const vector<int> &a, const vector<int> &b) {
    int need = a.size() + b.size() - 1;
    int nbase = 32 - __builtin_clz(need) - (need - need & (-need) == 0);
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
  vector<int> multiply_mod(const vector<int> &a, const vector<int> &b, int m, int eq = 0) {
    int need = a.size() + b.size() - 1;
    int nbase = 32 - __builtin_clz(need) - (need - need & (-need) == 0);
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
  vector<int> square_mod(const vector<int> &a, int m) {
    return multiply_mod(a, a, m, 1);
  }
}
```

### ntt

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
      if (rev) e = qpow(e, mod - 2, mod); // exgcd is better
      for (int i = 0; i < n; i += 2 * k) {
        ll w = 1;
        for (int j = 0; j < k; j++, w = w * e % mod) {
          int x = a[i + j], y = w * a[i + j + k] % mod;
          a[i + j] = (x + y) % mod, a[i + j + k] = (x - y + mod) % mod;
        }
      }
    }
    if (rev) {
      int inv = qpow(n, mod - 2, mod); // exgcd is better
      for (int i = 0; i < n; i++) a[i] = 1ll * a[i] * inv % mod;
    }
  }
  // mod = 998244353 = (119 << 23) + 1, root = 3, // = (119 << 23, 3)
  // For p < 2^30, (5 << 25, 3), (7 << 26, 3),
  // (479 << 21, 3) and (483 << 21, 5), last two are > 10^9
  vector<int> conv(const vector<int>& a, const vector<int>& b, const int mod = (119 << 23) + 1, int root = 3) {
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
```

<!-- 莫格马利取模 什么来的啊 有啥用啊= = -->

<!-- linear_seq 不会用啊 暂时不放 -->

<!-- 母函数模板 lhy说没用啊 没用就不放了啊 -->

