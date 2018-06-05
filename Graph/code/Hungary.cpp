#include <algorithm>
#include <vector>
#define szz(a) (int)a.size()
#define all(a) a.begin(),a.end()
using namespace std;

struct MaxMatch {
	int N, Ans;
	vector<vector<int>> e;
	vector<int> Link, Vis;

	MaxMatch(int n = 0) : N(n), e(n + 5), Link(n + 5), Vis(n + 5) { init(); }
	
	void init(const int& n = -1)
	{
		if (n != -1) N = n;
		for (int i = 0; i <= N; i++) e[i].clear();
	}

	inline void addedge(const int& u, const int& v)
	{
		e[u].push_back(v);
	}

	bool Find(const int& u)
	{
		for (int i = 0; i < szz(e[u]); i++) {
			int v = e[u][i];
			if (Vis[v]) continue;
			Vis[v] = 1;
			if (Link[v] == -1 || Find(Link[v])) {
				Link[v] = u;
				// Link[u] = v;
				return true;
			}
		}
		return false;
	}

	void Solve()
	{
		Ans = 0;
		fill(all(Link), -1);
		for (int i = 1; i <= N; i++) {
			if (Link[i] != -1) continue;
			fill(all(Vis), 0);
			if (Find(i)) ++Ans;
		}
	}
};
