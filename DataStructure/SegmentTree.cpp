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
		else if (mid > l) { 
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
		else if (mid > l) {
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