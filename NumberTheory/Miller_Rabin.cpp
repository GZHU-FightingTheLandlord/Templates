typedef long long ll;

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

ll qpow(ll a, ll t, ll mod) {
  ll b = 1;
  for (; t > 0; t >>= 1, a = a * a % mod) {
    if (t & 1) {
      b = b * a % mod;
    }
  }
  return b;
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

// WA -> make `t` higher
// TLE -> make `t` lower
bool isprime(ll n, int t=6) {
  if (n < 2) return false;
  if (n < psize) return !isp[n];
  if (~n & 1) return false;
  for (int j = 0; j < t; j++) {
    if (witness(rand() % (n - 1) + 1, n)) {
      return false;
    }
  }
  return true;
}
