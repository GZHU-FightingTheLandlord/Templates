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
        if (i == k) { y = x0, k <<= 1; }
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