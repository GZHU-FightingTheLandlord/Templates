#include <algorithm>
#include <string.h>
#include <vector>
using namespace std;

struct Dsu {
	int N, cnt;
	vector<int> Root;

	// N -> MaxNum
	Dsu(int num = 0) : N(num), Root(num + 5), cnt(num) {
		for (int i = 0; i < num; i++) Root[i] = i;
	}

	void init(int n = -1) // Initial
	{
		if (n != -1) N = n;
		cnt = N;
		for (int i = 0; i <= N; i++) Root[i] = i;
	}

	int find(int x) // Get Ancestor
	{
		if (x == Root[x]) return x;
		else return Root[x] = find(Root[x]);
	}

	bool Union(int a, int b)
	{
		if (!same(a, b)) {
			cnt--; // Trees in forest
			Root[Root[b]] = Root[a]; // Link b to a
			return true;
		}
		return false;
	}

	bool same(int a, int b)
	{
		a = find(a);
		b = find(b);
		return a == b; // Same Ancestor?
	}
};

struct Kruskal {

#define all(a) a.begin(),a.end()
	
	struct edge{
		int u, v, w;
		edge(int uu = 0, int vv = 0, int ww = 0) {
			u = uu, v = vv, w = ww;
		}
		bool operator< (const edge& b)const {
			return w < b.w;
		}
	};

	int N, ans, ecnt;
	vector<edge> e;
	vector<int> used;
	Dsu D;

	Kruskal(int MAX = 0) : N(MAX), e(MAX * MAX / 2 + 5), D(MAX), used(MAX + 5) {}
	
	void init(int n = -1)
	{
		if (n != -1) N = n;
		ecnt = 0; D.init(N);
	}

	void addedge(int u, int v, int w)
	{
		e[ecnt++] = edge(u, v, w);
	}

	void solve()
	{
		ans = 0;
		sort(all(e));
		fill(all(used), 0);
		for (int i = 0; i < ecnt; i++) {
			if (!D.same(e[i].u, e[i].v)) {
				used[i] = 1;
				ans += e[i].w;
				D.Union(e[i].v, e[i].u);
			}
		}
	}
};