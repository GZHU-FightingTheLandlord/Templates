#include <algorithm>
#include <queue>
#include <vector>
using namespace std;

template <typename T> struct Dinic {
	int N, cnt;
	int st, ed;
	T Ans, temp;

	static const ll inf = 0x3f3f3f3f3f3f3f3f;

	struct edge {
		int v, next;
		T w;
		edge(const int& vv = 0, const T& ww = 0, const int& nn = -1) {
			v = vv, w = ww, next = nn;
		}
	};

	vector<edge> e;
	vector<int> head, level, cur;

	Dinic(const int& n = 0, const int& eg = 0, const int& s = 1, const int& t = 0x3f3f3f3f) : e(eg + 5), head(n + 5), level(n + 5) {
		N = n, st = s, ed = min(n, t);
		init();
	}

	inline void init()
	{
		cnt = 0;
		fill(head.begin(), head.end(), -1);
	}

	inline void addedge(const int& u, const int& v, const T& w)
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
		return level[ed];
	}

	T Dfs(T flow, const int& u)
	{
		if (u == ed || flow == 0) return flow;

		T delta = 0;
		//for (int& i = cur[u]; i != -1 && flow > 0; i = e[i].next) {
		for (int i = head[u]; i != -1 && flow > 0; i = e[i].next) {
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
		if (delta == 0) level[u] = -1;
		return delta;
	}

	T solve()
	{
		Ans = 0;
		while (Bfs()) {
			//cur = head;
			while ((temp = Dfs(inf, st)) > 0) {
				Ans += temp;
			}
		}
		return Ans;
	}
};
