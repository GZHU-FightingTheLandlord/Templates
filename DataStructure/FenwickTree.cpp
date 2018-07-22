#include <stdio.h>
#include <algorithm>
#include <vector>
using namespace std;

/*
	New FenwickTree Solution!!
	BIT::Size -> MAXIMUM
	init -> reset [0, len]
	upd -> val[i] += dlt
	gsum -> sum_z [0, i]
	query -> sum_z [l, r]
*/
template <typename T> struct BIT {
	int Size;
	vector<T> base;
	BIT(int n = 0) : Size(n), base(n + 5, 0) {}
#define lowbit(i) (~i & i + 1)
	void init(int len = 0) { for (; len >= 0; len--) base[len] = 0; }
	void upd(int i, T dlt) { for (; i <= Size; i += lowbit(i)) base[i] += dlt; }
	T gsum(int i) { T ret = 0; for (; i >= 0; i -= lowbit(i)) ret += base[i]; return ret; }
	T query(int l, int r) { T ret = 0; for (; r >= 0; r -= lowbit(r)) ret += base[r]; for (--l; l >= 0; l -= lowbit(l)) ret -= base[l]; return ret; }
};

template <typename T> struct Fenwick {
	int n;
	vector<T> sum_;

#define lowbit(x) (x & (-x))

	Fenwick(int n_ = 0) : n(n_), sum_(n + 5, 0) {}

	void init(int n_ = -1)
	{
		if (n_ != -1) n = n_;
		for (int i = 0; i <= n; i++) sum_[i] = 0;
	}

	void upd(int i, T delta)
	{
		for (; i <= n; i += lowbit(i)) sum_[i] += delta;
	}

	T getsum(int i)
	{
		T res = 0;
		for (; i > 0; i -= lowbit(i)) res += sum_[i];
		return res;
	}

	T query(int l, int r)
	{
		return getsum(r) - getsum(l - 1);
	}
};
