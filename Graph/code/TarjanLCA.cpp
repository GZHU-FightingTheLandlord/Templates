#include <algorithm>
#include <vector>
using namespace std;

struct Dsu {
	int N, cnt;
	vector<int> Root;

	Dsu(int num = 0) : N(num), Root(num + 5), cnt(num) {
		for (int i = 0; i < num; i++) Root[i] = i;
	}

	void init() // Initial
	{
		cnt = N;
		for (int i = 0; i < N; i++) Root[i] = i;
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

struct TarjanLCA {
	int N;
	Dsu D;
	vector<vector<int>> e;
	vector<int> Color;

	struct Node {
		int u, v, w;
		Node(int uu = 0, int vv = 0, int ww = 0) {
			u = uu, v = vv, w = ww;
		}
	};

	vector<vector<Node>> Q;
	vector<Node> Ans;

	TarjanLCA(int n = 0) : N(n + 5), D(n + 5), \
		e(n + 5), Color(n + 5), Q(n + 5) { init(); }

	void init(int n = -1) // Initial
	{
		if (n < 0) n = N - 1;
		for (int i = 1; i <= n; i++) e[i].clear(), Q[i].clear();
		fill(Color.begin(), Color.end(), 0);
		Ans.clear();
		D.init();
	}

	void addedge(const int& u, const int& v)
	{
		e[u].push_back(v);
		e[v].push_back(u);
	}

	void addquery(const int& u, const int& v)
	{
		Q[u].emplace_back(u, v, Ans.size());
		Q[v].emplace_back(v, u, Ans.size());
		Ans.emplace_back(u, v, -1);
	}

	void Dfs(const int& u, const int& pre)
	{
		for (int i = 0; i < e[u].size(); i++) {
			if (Color[e[u][i]] || e[u][i] == pre) continue;
			Dfs(e[u][i], u); // Dfs child
		}

		for (int i = 0; i < Q[u].size(); i++) {
			if (!Color[Q[u][i].v]) continue;
			Ans[Q[u][i].w].w = D.find(Q[u][i].v); // Get LCA
		}

		Color[u] = 1;
		if (pre != -1) D.Union(pre, u); // Attention!! Link u to pre!!!
	}

	void Solve()
	{
		Dfs(1, -1); // Consider 1 as the highest ancestor
	}
};
