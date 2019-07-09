ll fpow(ll a, ll b, ll c) {
  ll ans = 1;
  for(a %= c; b; a = a * a % c, b >>= 1) {
    if(b & 1) {
      ans = ans * a % c;
    }
  }
  return ans;
}

const int MOD = 1e9 + 7;
const int inv2 = fpow(2, MOD - 2, MOD);

int trans(int n) {
  return 1 << (32 - __builtin_clz(n) - (n - (n & (-n)) == 0));
}

void fwt(vector<int> &a) {
  const int n = trans(a.size());
  a.resize(n);
  for (int d = 1; d < n; d <<= 1) {
    for (int i = 0, k = d << 1; i < n; i += k) {
      for (int j = 0; j < d; j++) {
        int x = a[i + j], y = a[i + j + d];
        a[i + j] = (x + y) % MOD, a[i + j + d] = (x - y + MOD) % MOD;
        // xor : a[i + j] = x + y, a[i + j + d] = x - y
        // and : a[i + j] = x + y
        // or : a[i + j + d] = x + y
      }
    }
  }
}

void ifwt(vector<int> &a) {
  const int n = trans(a.size());
  a.resize(n);
  for (int d = 1; d < n; d <<= 1) {
    for (int i = 0, k = d << 1; i < n; i += k) {
      for (int j = 0; j < d; j++) {
        int x = a[i + j], y = a[i + j + d];
        a[i + j] = 1ll * (x + y) * inv2 % MOD, a[i + j + d] = (1ll * (x - y) * inv2 % MOD + MOD) % MOD;
        // xor : a[i + j] = (x + y) / 2, a[i + j + d] = (x - y) / 2
        // and : a[i + j] = x - y
        // or : a[i + j + d] = y - x;
      }
    }
  }
}
