pll crt(const vector<pll> & v) {
    //v里每个pll中first为被模数，second为模数
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