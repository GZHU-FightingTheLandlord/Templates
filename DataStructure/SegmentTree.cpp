#include <algorithm>
#include <vector>
using namespace std;

template <typename T> struct SegmentTree {
	struct Node	{
		int l, r;
		T val, change, cover;
		Node() {}
	};

	int N;
	vector<Node> v;
	vector<bool> covered;

	SegmentTree(int n) : N(n), v(n << 2), covered(n << 2) {}

	void build(int index, int l, int r, const T arr[])
	{
		v[index].l = l, v[index].r = r, v[index].change = 0, v[index].cover = 0;
		if (l == r) {
			v[index].val = arr[l];
			return;
		}
		int mid = (l + r) >> 1;
		build(index << 1, l, mid, arr);
		build(index << 1 | 1, mid + 1, r, arr);
		v[index].val = min(v[index << 1].val, v[index << 1 | 1].val);
	}

	void pushdown(const int &index)
	{
		if (covered[index])	{
			v[index].val = v[index].cover;
			v[index << 1].cover = v[index << 1 | 1].cover = v[index].cover;
			covered[index << 1] = covered[index << 1 | 1] = 1;
			covered[index] = 0;
		}
		if (v[index].change) {
			v[index].val += v[index].change;
			v[index << 1].change += v[index].change;
			v[index << 1 | 1].change += v[index].change;
			v[index].change = 0;
		}
	}

	void update(int index, int l, int r, const T &delta)
	{
		if (l <= v[index].l && v[index].r <= r) {
			v[index].change += delta;
			return;
		}
		pushdown(index);

		int mid = (v[index].l + v[index].r) >> 1;
		if (r <= mid) {
			update(index << 1, l, r, delta);
		}
		else if (mid < l) { 
			update(index << 1 | 1, l, r, delta);
		}
		else {
			update(index << 1, l, r, delta);
			update(index << 1 | 1, l, r, delta);
			v[index].val = min(query(index << 1, v[index].l, mid), query(index << 1 | 1, mid + 1, v[index].r));
		}
	}

	void modify(int index, int l, int r, const T &delta)
	{
		if (l <= v[index].l && v[index].r <= r) {
			v[index].change = 0;
			v[index].cover = delta;
			covered[index] = 1;
			return;
		}
		pushdown(index);

		int mid = (v[index].l + v[index].r) >> 1;
		if (r <= mid) {
			modify(index << 1, l, r, delta);
		}
		else if (mid < l) {
			modify(index << 1 | 1, l, r, delta);
		}
		else {
			modify(index << 1, l, r, delta);
			modify(index << 1 | 1, l, r, delta);
			v[index].val = min(query(index << 1, v[index].l, mid), query(index << 1 | 1, mid + 1, v[index].r));
		}
	}

	T query(int index, int l, int r)
	{
		if (v[index].l == l && v[index].r == r) {
			if (covered[index])	{
				return v[index].cover;
			}
			return v[index].val + v[index].change;
		}
		pushdown(index);
		int mid = (v[index].l + v[index].r) >> 1;
		if (r <= mid) {
			return query(index << 1, l, r);
		}
		else if (mid < l) {
			return query(index << 1 | 1, l, r);
		}
		else {
			return min(query(index << 1, l, mid), query(index << 1 | 1, mid + 1, r));
		}
	}
};

// ********************************************

const int MAX = 1e5 + 5;
struct node {
	int l, r, val;
	int lazy;
}v[MAX << 2];

inline void pushup(int i) {
	v[i].val = min(v[i << 1].val, v[i << 1 | 1].val);
}

inline void pushdown(int i) {
	if (v[i].lazy) {
		int lson = i << 1, rson = i << 1 | 1;
		v[lson].val += v[i].lazy; v[rson].val += v[i].lazy;
		v[lson].lazy += v[i].lazy; v[rson].lazy += v[i].lazy;
		v[i].lazy = 0;
	}
}

void build(int *arr, int i, int l, int r) {
	v[i].l = l, v[i].r = r, v[i].lazy = 0;
	if (l == r) {
		v[i].val = arr[l];
		return;
	}
	int mid = (l + r) / 2;
	build(arr, i << 1, l, mid); build(arr, i << 1 | 1, mid + 1, r);
	pushup(i);
}

void update(int i, int L, int R, int k) {
	if (L <= v[i].l && v[i].r <= R) {
		v[i].val += k, v[i].lazy += k;
		return;
	}
	pushdown(i);
	int mid = (v[i].l + v[i].r) / 2;
	if (R <= mid) update(i << 1, L, R, k);
	else if (L > mid) update(i << 1 | 1, L, R, k);
	else update(i << 1, L, R, k), update(i << 1 | 1, L, R, k);
	pushup(i);
}

int query(int i, int L, int R) {
	if (L <= v[i].l && v[i].r <= R) {
		return v[i].val;
	}
	pushdown(i);
	int mid = (v[i].l + v[i].r) / 2;
	if (R <= mid) return query(i << 1, L, R);
	else if (L > mid) return query(i << 1 | 1, L, R);
	else return min(query(i << 1, L, R), query(i << 1 | 1, L, R));
}
