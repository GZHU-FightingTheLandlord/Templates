const int maxn = 1e5 + 5;

int a[maxn];

// lowbit(i) = i & -i
// lowbit(i) = ~i & i + 1

void upd(int i, int x) {
	for (; i < maxn; i += (i & (-i))) {
		a[i] += x;
	}
}

int query(int i) {
	int res = 0;
	for (; i > 0; i -= (i & (-i))) {
		res += a[i];
	}
	return res;
}

inline int query(int l, int r) {
	return query(r) - query(l - 1);
}

// 以x为上界，找最后一个满足\sum_{i = 1}^{k}的点k
int upper_bound(int x) {
    int res = 0, ptr = 31 - __builtin_clz(n);
    for (; ~ptr; --ptr) {
        int p = res | (1 << ptr);
        if (p <= n && a[p] <= x) {
            x -= a[p];
            res |= 1 << ptr;
        }
    }
    return res;
}