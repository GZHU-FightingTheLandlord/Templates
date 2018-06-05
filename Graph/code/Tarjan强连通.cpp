#include <algorithm>
#include <string.h>
#include <vector>
#include <stack>
#define szz(a)	(int)a.size()
using namespace std;

struct Tarjan {
	int N, idx, cnt;
	vector<vector<int>> e;
	stack<int> S;
	vector<int> def, low, ins, sccno;

	Tarjan(int num = 0) : N(num), e(num + 5), def(num + 5), \
		low(num + 5), ins(num + 5, 0), sccno(num + 5) { init(); }

	void init(int n = -1)
	{
		if (n != -1) N = n;
		for (int i = 0; i <= N; i++) e[i].clear(), sccno[i] = -1;
	}

	inline void addedge(const int& u, const int& v)
	{
		e[u].push_back(v);
	}

	void Dfs(int u)
	{
		def[u] = ++idx;
		low[u] = idx;
		S.push(u); ins[u] = 1;

		int v;
		for (int i = 0; i < szz(e[u]); i++) {
			v = e[u][i];
			if (!def[v]) {
				Dfs(v);
				low[u] = min(low[u], low[v]);
			}
			else if (ins[v]) {
				low[u] = min(low[u], def[v]);
			}
		}

		if (def[u] == low[u]) {
			++cnt;
			do {
				v = S.top();
				S.pop(); ins[v] = 0;
				sccno[v] = cnt;
			} while (u != v);
		}
	}

	void solve()
	{
		idx = 0; cnt = 0;
		fill(def.begin(), def.end(), 0);
		fill(low.begin(), low.end(), 0);
		while (!S.empty()) S.pop();

		for (int i = 1; i <= N; i++) if (!def[i]) Dfs(i);
	}
};
