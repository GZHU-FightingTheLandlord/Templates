#include <algorithm>
#include <queue>
#include <vector>
using namespace std;

template <typename T> struct Dinic {
	int N, cnt;
	int st, ed;
	T Ans;

	static const int inf = 0x3f3f3f3f;

	struct edge {
		int v, next;
		T w;
		edge(int vv = 0, T ww = 0, int nn = -1) {
			v = vv, w = ww, next = nn;
		}
	};

	vector<edge> e;
	vector<int> head, level, cur;

	Dinic(int n = 0, int eg = 0, int s = 1, int t = 0x3f3f3f3f) : e(eg + 5), head(n + 5), level(n + 5) {
		N = n, st = s, ed = min(n, t);
		init();
	}

	void init()
	{
		cnt = 0;
		fill(head.begin(), head.end(), -1);
	}

	void addedge(int u, int v, int w) 
	{
		e[cnt] = edge(v, w, head[u]); head[u] = cnt++;
		e[cnt] = edge(u, 0, head[v]); head[v] = cnt++;
	}

	bool Bfs()
	{
		fill(level.begin(), level.end(), 0);
		queue<int> Que;

		Que.push(st); level[st] = 1;
		while (!Que.empty()) {
			int u = Que.front();
			Que.pop();

			for (int i = head[u]; i != -1; i = e[i].next) {
				int v = e[i].v;
				if (!level[v] && e[i].w > 0) {
					level[v] = level[u] + 1;
					Que.push(v);
				}
			}
		}
		
		return level[ed] != 0;
	}

	T Dfs(T flow, int u)
	{
		if (u == ed || flow == 0) return flow;

		T delta = 0;
		for (int& i = cur[u]; i != -1 && flow > 0; i = e[i].next) {
			int v = e[i].v;
			if (level[v] == level[u] + 1 && e[i].w > 0) {
				T tmp = Dfs(min(e[i].w, flow), v);
				if (tmp == 0) continue;
				e[i].w -= tmp;
				e[i ^ 1].w += tmp;
				delta += tmp;
				flow -= tmp;
			}
		}

		return delta;
	}

	T solve()
	{
		Ans = 0;
		while (Bfs()) {
			cur = head;
			Ans += Dfs(inf, st);
		}
		return Ans;
	}
};
