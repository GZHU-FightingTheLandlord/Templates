pll exgcd(const long long x, const long long y) {
    if (!y) {
        return make_pair(1, 0);
    }
    pll cur = exgcd(y, x % y);
    return make_pair(cur.second, cur.first - (x / y) * cur.second);
}
