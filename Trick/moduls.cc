typedef long long ll;
const int MOD = 998244353;

inline ll add(ll a, ll b) {
    a += b;
    return (a < MOD) ? a : (a - MOD);
}

inline ll sub(ll a, ll b) {
    a -= b;
    return (a < 0) ? (a + MOD) : a;
}

inline ll mul(ll a, ll b) {
    a *= b;
    return a < MOD ? a : (a - a / MOD * MOD);
}