typedef long long ll;

// mod <= 1e12
inline ll mul(ll a, ll b, ll mod) {
    return (((a * (b >> 20) % mod) << 20) + (a * (b & ((1 << 20) - 1)))) % mod;
}

inline ll mul(ll a, ll b, ll mod) {
    ll d = (ll)floor(a * (long double)b / mod + 0.5);
    ll ret = (a * b - d * mod) % mod;
    if (ret < 0) ret += mod;
    return ret;
}
