// C(n, m) % p = C(n % p, m % p) * C(n / p, m / p) % p
// p is a prime

typedef long long ll;
const int MOD = 1e9 + 7;

// calc C(n, m) % MOD
inline ll C(ll n, ll m);

ll lucas(ll n, ll m) {
    if (n < MOD && m < MOD) {
        return C(n, m);
    }
    return C(n % MOD, m % MOD) * lucas(n / MOD, m / MOD) % MOD;
}