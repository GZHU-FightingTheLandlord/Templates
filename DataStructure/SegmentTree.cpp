#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

const int MAXN = 1e5 + 5;

#define lson(x) ((x) << 1)
#define rson(x) ((x) << 1 | 1)

struct SegNode{
	int l, r;
	ll val;
	ll lazy;	// 延迟标记， 表示区域内每个单独节点的增量
	ll force;	// 覆盖标记， 表示区域内每个单独节点被修改后的值
}node[MAXN << 2]; // root为1， 四倍空间初始
bool covered[MAXN << 2]; // 节点区域是否被覆盖

// build(1, 1, n, arr)
// arr为初值， 以1为root， 建立区间为[1, n]的线段树
void build(int idx, int l, int r, ll *arr)
{
	node[idx].l = l, node[idx].r = r;
	node[idx].lazy = node[idx].force = 0;
	if (l == r)
	{
		node[idx].val = arr[l];
	}
	else
	{
		int mid = (l + r) / 2;
		build(lson(idx), l, mid, arr);
		build(rson(idx), mid + 1, r, arr);
		node[idx].val = node[lson(idx)].val + node[rson(idx)].val;
	}
}

// 下传标记
void pushdown(int idx)
{
	if (covered[idx])
	{
		node[idx].val = node[idx].force * (node[idx].r - node[idx].l + 1);
		covered[lson(idx)] = covered[rson(idx)] = true;

		node[lson(idx)].lazy = node[rson(idx)].lazy = 0;
		node[lson(idx)].force = node[rson(idx)].force = node[idx].force;
		covered[idx] = false;
	}
	if (node[idx].lazy)
	{
		node[idx].val += node[idx].lazy * (node[idx].r - node[idx].l + 1);

		node[lson(idx)].lazy += node[idx].lazy;
		node[rson(idx)].lazy += node[idx].lazy;
		node[idx].lazy = 0;
	}
}


// query(1, l, r)
// 查询区间[l, r]的sum
ll query(int idx, int l, int r)
{
	if (node[idx].l == l && node[idx].r == r)
	{
		if (covered[idx])
			return (node[idx].force + node[idx].lazy) * (node[idx].r - node[idx].l + 1);
		return node[idx].val + node[idx].lazy * (node[idx].r - node[idx].l + 1);
	}
	pushdown(idx);

	int mid = (node[idx].l + node[idx].r) / 2;
	if (r <= mid)
		return query(lson(idx), l, r);
	if (mid < l)
		return query(rson(idx), l, r);
	return query(lson(idx), l, mid) + query(rson(idx), mid + 1, r);
}


// update(1, l, r, diff)
// 将[l, r]区间内元素增加diff
void update(int idx, int l, int r, int diff)
{
	if (l <= node[idx].l && node[idx].r <= r)
	{
		node[idx].lazy += diff;
		return ;
	}
	pushdown(idx);

	int mid = (node[idx].l + node[idx].r) / 2;
	if (r <= mid)
		update(lson(idx), l, r, diff);
	else if (mid < l)
		update(rson(idx), l, r, diff);
	else
	{
		update(lson(idx), l, r, diff);
		update(rson(idx), l, r, diff);
	}
	node[idx].val = query(lson(idx), node[idx].l, mid) + query(rson(idx), mid + 1, node[idx].r); 
}

// modify(1, l, r, val);
// 将[l, r]内的元素修改为val
void modify(int idx, int l, int r, int val)
{
	if (l <= node[idx].l && node[idx].r <= r)
	{
		node[idx].lazy = 0;
		node[idx].force = val;
		covered[idx] = true;
		return ;
	}
	pushdown(idx);

	int mid = (node[idx].l + node[idx].r) / 2;
	if (r <= mid)
		modify(lson(idx), l, r, val);
	else if (mid < l)
		modify(rson(idx), l, r, val);
	else
	{
		modify(lson(idx), l, r, val);
		modify(rson(idx), l, r, val);
	}
	node[idx].val = query(lson(idx), node[idx].l, mid) + query(rson(idx), mid + 1, node[idx].r);
}


