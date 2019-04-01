// https://www.spoj.com/problems/MOD/en/
// https://vjudge.net/problem/SPOJ-MOD
// 返回最小正整数n使得 a^n mod m = b; O(sqrt(m))

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
