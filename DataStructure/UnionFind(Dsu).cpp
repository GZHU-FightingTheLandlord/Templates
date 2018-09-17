#include <bits/stdc++.h>
using namespace std;

struct Dsu {
	static const int maxn = 1e5 + 5;
	int fa[maxn], sz[maxn];
	void init(int n) {
		for (int i = 0; i <= n; i++) {
			fa[i] = i, sz[i] = 0;
		}
	}
	int find(int x) {
		return x == fa[x] ? x : fa[x] = find(fa[x]);
	}
	bool unite(int u, int v) {
		int a = find(u), b = find(v);
		if (a == b) return false;
		if (sz[a] < sz[b]) fa[a] = b;
		else fa[b] = a, sz[a] += (sz[a] == sz[b]);
		return true;
	}
};