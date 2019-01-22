// https://www.spoj.com/problems/MOD/en/
// https://vjudge.net/problem/SPOJ-MOD
// a^x mod p = b; find minimum x; O(sqrt(p));
// 可能会卡常的地方： unordered_map, fmul


#include <cstdio>
#include <cmath>
#include <unordered_map>
#include <cstring>
using namespace std;

inline long long fmul(long long a, long long c, const long long &Mod = 998244353ll) {
    a %= Mod, c %= Mod;
    return (a * c - (long long)(((long double)a * c + 0.5) / Mod) * Mod) % Mod;
}

long long gcd(const long long x, const long long y) {
    return !y ? x : gcd(y, x % y);
}

long long fpow(long long a, long long b, const long long &p = 998244353ll) {
    long long ans = 1ll;
    for (a %= p; b; a = fmul(a, a, p), b >>= 1) {
        if (b & 1) {
            ans = fmul(ans, a, p);
        }
    }
    return ans;
}

// a^x = b (mod p), with a, b, p return x
long long exbsgs(long long a, long long b, long long p) {
    a %= p, b %= p;
    if (a == 0) {
        long long ret = b > 1 ? -1 : (long long)(b == 0 && p > 1);
        if(ret == -1) {
            return ret;
        }
    }
    long long g, c = 0, q = 1;
    while ((g = gcd(a, p)) != 1) {
        if (b == 1) {
            return c;
        }
        if (b % g) {
            return -1;
        }
        ++c;
        b /= g, p /= g;
        q = fmul(a / g, q, p);
    }
    unordered_map<long long, long long> mp;
    long long m = sqrt(p);
    for (long long i = 1, t = fmul(b, a, p); i <= m; i++, t = fmul(t, a, p)) {
        mp[t] = i;
    }
    for (long long i = m, t = fpow(a, m, p); i - m < p - 1; i += m) {
        q = fmul(q, t, p);
        if (mp.count(q)) {
            return i - mp[q] + c;
        }
    }
    return -1;
}

int main() {
    long long a, b, c;
    while(~scanf("%lld %lld %lld", &a, &b, &c) && a + b + c) {
        long long ans = exbsgs(a, c, b);
        if(ans == -1) {
            printf("No Solution\n");
        } else {
            printf("%lld\n", ans);
        }
    }
}
