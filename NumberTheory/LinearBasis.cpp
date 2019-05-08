struct LinearBasis {
    const static int MAXL = 50;
    long long a[MAXL + 1];
    LinearBasis() {
        memset(a, 0, sizeof a);
    }
    void insert(long long t) {
        for (int j = MAXL; j >= 0; j--) {
            if (!(t & (1ll << j))) continue;
            if (a[j]) t ^= a[j];
            else {
                for (int k = 0; k < j; k++) if (t & (1ll << k)) t ^= a[k];
                for (int k = j + 1; k <= MAXL; k++) if (a[k] & (1ll << j)) a[k] ^= t;
                a[j] = t;
                return;
            }
        }
    }
};
