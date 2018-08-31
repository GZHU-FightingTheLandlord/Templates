#include <algorithm>
#include <string.h>
#include <vector>
using namespace std;

namespace Factor {
    using pll = pair<int64_t, int64_t>;

    const int maxn = 1e6 + 5;

    bool isp[maxn], ok_init;
    int prime[maxn / 10], tot_prime, phi[maxn];

    inline void init() {
        register int i, j;
        ok_init = true;
        memset(isp, true, sizeof isp);
        for (i = 2, isp[0] = isp[1] = false, phi[1] = 1; i < maxn; i++) {
            if (isp[i]) {
                prime[tot_prime++] = i, phi[i] = i - 1;
            }
            for (j = 0; j < tot_prime && i * prime[j] < maxn; j++) {
                isp[i * prime[j]] = false;
                if (i % prime[j]) {
                    phi[i * prime[j]] = phi[i] * (prime[j] - 1);
                }
                else {
                    phi[i * prime[j]] = phi[i] * prime[j];
                    break;
                }
            }
        }
    }

    int64_t mul(int64_t a, int64_t b, int64_t p) {
        if (p <= (int)1e9) return a * b % p;
        else if (p <= 1000000000000ll) {
            return (((a * (b >> 20) % p) << 20) + (a * (b & ((1 << 20) - 1)))) % p;
        }
        else {
            int64_t d = (int64_t)floor(a * (long double)b / p + 0.5);
            int64_t ret = (a * b - d * p) % p;
            if (ret < 0) ret += p;
            return ret;
        }
    }

    int64_t qpow(int64_t a, int t, int64_t p) {
        int64_t b = 1;
        for (; t > 0; t >>= 1, a = mul(a, a, p)) {
            if (t & 1) {
                b = mul(b, a, p);
            }
        }
        return b;
    }

    inline int* prime_table() {
        if (!ok_init) init();
        return prime;
    }
    
    inline int* phi_table() {
        if (!ok_init) init();
        return phi;
    }

    int64_t eular(int64_t x) {
        if (!ok_init) init();
        if (x < maxn) return phi[x];
        int64_t i, res = x;
        for (i = 0; i < tot_prime && prime[i] * prime[i] <= x; ++i) {
            if (x % prime[i] == 0) {
                res = res / i * (i - 1);
                while (x % prime[i] == 0) x /= prime[i];
            }
        }
        return (x > 1) ? (res / x * (x - 1)) : res;
    }

    vector<pll> prifac;
    vector<int64_t> facres;

    inline void deal() {
        sort(prifac.begin(), prifac.end());
        int i, j;
        for (i = 1, j = 0; i < prifac.size(); i++) {
            if (prifac[i].first == prifac[j].first) {
                prifac[j].second += prifac[j].second;
            }
            else prifac[++j] = prifac[i];
        }
        prifac.resize(j + 1);
    }

    void DFS(int64_t x, int i) {
        if (i == prifac.size()) {
            facres.push_back(x);
            return;
        }
        DFS(x, i + 1);
        for (int j = 1; j <= prifac[i].second; j++) {
            DFS(x *= prifac[i].first, i + 1);
        }
    }

    vector<pll> prifactor(int64_t x) {
        if (!ok_init) init();
        prifac.clear();
        for (int i = 0; i < tot_prime && 1ll * prime[i] * prime[i] <= x; i++) {
            if (x % prime[i] == 0) {
                prifac.emplace_back(prime[i], 0);
                while (x % prime[i] == 0) {
                    x /= prime[i];
                    prifac.back().second++;
                }
            }
        }
        if (x > 1) prifac.emplace_back(x, 1);
        deal();
        return prifac;
    }

    vector<int64_t> factor(int64_t x) {
        prifactor(x);
        facres.clear();
        DFS(1, 0);
        sort(facres.begin(), facres.end());
        return facres;
    }
};