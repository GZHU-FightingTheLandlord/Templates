typedef long long ll;

int tot;
ll factor[10000];

ll pollard_rho(ll x, ll c) {
  ll i = 1, k = 2;
  ll x0 = rand() % (x - 1) + 1, y = x0;
  while (true) {
    i++;
    x0 = (mul(x0, x0, x) + c) % x;
    ll d = __gcd(y - x0, x);
    if(d < 0) d = -d;
    if (d != 1 && d != x) return d;
    if (y == x0) return x;
    if (i == k) y = x0, k <<= 1;
  }
}

void findfac(ll n, int k=107) {
  if(n == 1) return;
  if (isprime(n)) {
    factor[tot++] = n;
    return;
  }
  ll p = n;
  int c = k;
  while (p >= n) {
    p = pollard_rho(p, c--);
  }
  findfac(p, k), findfac(n / p, k);
}
