const int maxn = 1e5 + 5;

int a[maxn];

inline void init(int n) {
	for (; ~n; --n) a[n] = 0;
}

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

int upper_bound(int x) {
	int res = 0, ptr = 0;
	while ((1 << ptr + 1) <= n) ++ptr;
	for (; ~ptr; --ptr) {
		int p = res + (1 << ptr);
		if (p <= n && a[p] <= x) {
			x -= a[p];
			res += 1 << ptr;
		}
	}
	return res;
}