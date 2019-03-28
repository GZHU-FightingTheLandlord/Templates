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

// https://codeforces.com/gym/100633/problem/J

/* description:
    typedef long long ll;
    typedef pair<long long, long long>pll;
    mod(x, p) return ((x % p) + p) % p
    fpow(a, b, p) return a^b mod p
    exgcd(a, b) return x, y that make x * a + y * b = gcd(x, y)
*/
