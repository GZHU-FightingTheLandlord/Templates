#include <bits/stdc++.h>
using namespace std;

/*
    it::upper_bound(x)
    O(log n) solution
    find the largest prefix with sum not exceeding x.
*/

struct it {
    const int it_max = 1 << 18;
    int n, a[it_max];

    void init(int k = 0) {
        for (n = k; k; a[k--] = 0);
    }
    
    #define lowbit(x) (x & (-x))

    void upd(int i, int x) {
        for (; i <= n; i += lowbit(i)) {
            a[i] += x;
        }
    }
    int query(int l, int r) {
        int ret = 0;
        for (; r > 0; r -= lowbit(r)) {
            ret += a[r];
        }
        for (--l; l > 0; l -= lowbit(l)) {
            ret -= a[l];
        }
        return ret;
    }
    int upper_bound(int x) {
        int res = 0, ptr = 0;
        while ((1 << (ptr + 1)) <= n) ++ptr;
        for (; ptr >= 0; ptr--) {
            int p = res + (1 << ptr);
            if (p <= n && a[p] <= x) {
                x -= a[p];
                res += 1 << ptr;
            }
        }
        return res;
    }
}bit;

int32_t main() {
    int n; cin >> n;
    bit.init(n);
    for (int i = 1; i <= n; i++) {
        int k; cin >> k;
        bit.upd(i, k);
    }
    int q; cin >> q;
    while (q--) {
        int t; cin >> t;
        if (t == 1) {
            // val[i] += x;
            int i, x; cin >> i >> x;
            bit.upd(i, x);
        }
        else if (t == 2) {
            // sum of [l, r]
            int l, r; cin >> l >> r;
            cout << bit.query(l, r) << endl;
        }
        else {
            // find the largest prefix with sum not exceeding x.
            int x; cin >> x;
            cout << bit.upper_bound(x) << endl;
        }
    }
}