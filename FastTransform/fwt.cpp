const int MOD = 1e9 + 7;

int qpow(int a, int t) {
    int b = 1;
    while (t > 1) {
        if (t & 1) b = b * a % MOD;
        a = a * a % MOD;
        t >>= 1;
    }
    return b;
}

void fwt(int a[], int n) {
    for (int d = 1; d < n; d <<= 1) {
        for (int i = 0, k = d << 1; i < n; i += k) {
            for (int j = 0; j < d; j++) {
                int x = a[i + j], y = a[i + j + d];
                a[i + j] = (x + y) % MOD;
                a[i + j + d] = (x - y + MOD) % MOD;
                // xor : a[i + j] = x + y, a[i + j + d] = x - y
                // and : a[i + j] = x + y
                // or : a[i + j + d] = x + y
            }
        }
    }
}

const int inv2 = qpow(2, MOD - 2);
void ifwt(int a[], int n) {
    for (int d = 1; d < n; d <<= 1) {
        for (int i = 0, k = d << 1; i < n; i += k) {
            for (int j = 0; j < d; j++) {
                int x = a[i + j], y = a[i + j + d];
                a[i + j] = 1ll * (x + y) * inv2 % MOD;
                a[i + j + d] = (1ll * (x - y) * inv2 % MOD + MOD) % MOD;
                // xor : a[i + j] = (x + y) / 2, a[i + j + d] = (x - y) / 2
                // and : a[i + j] = x - y
                // or : a[i + j + d] = y - x;
            }
        }
    }
}