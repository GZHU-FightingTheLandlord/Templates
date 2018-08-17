#include <algorithm>
#include <vector>
#include <queue>
using namespace std;

template <typename T> struct Node {
	int u;
	T w;
	Node(int u_ = 0, T w_ = 0) : u(u_), w(w_) {}
	// Dijkstra Priority_queue
	bool operator< (const Node& b) const {
		return w > b.w;
	}
};

template <typename T> struct Dijkstra {
	int N, st, ed;
	T Ans, inf;
	vector<vector<pair<int, T> > > e;
	vector<T> Dis;
	vector<bool> Vis;

	// num -> MaxNum
	Dijkstra(int num = 0) : N(num), e(num + 5), Dis(num + 5), Vis(num + 5) {
		Ans = 0; inf = 0x3f3f3f3f;
	}

	// 0 -> Don't reset the Graph
	// 1 -> Reset all
	void init(bool type = 1)
	{
		fill(Dis.begin(), Dis.end(), inf);
		fill(Vis.begin(), Vis.end(), 0);
		if (type) {
			e = vector<vector<pair<int, T> > >(N + 5);
		}
	}

	// from u to v, cost w
	void addedge(const int& u, const int& v, const T& w)
	{
		e[u].emplace_back(v, w);
		//e[v].emplace_back(u, w);
	}

	// Set the Starting Point st and Ending Point ed before using
	T solve()
	{
		init();
		priority_queue<Node<T> > Que;
		Que.emplace(st, 0); Dis[st] = 0;
		while (!Que.empty()) {
			int u = Que.top().u;
			Que.pop();

			if (ed != -1 && u == ed) return Ans = Dis[ed];
			if (Vis[u]) continue;
			Vis[u] = 1;

			for (int i = 0; i < e[u].size(); i++) {
				int v = e[u][i].first, w = e[u][i].second;
				if (!Vis[v] && Dis[v] > Dis[u] + w) {
					Dis[v] = Dis[u] + w;
					Que.emplace(v, Dis[v]);
				}
			}
		}
		return Ans = -1;
	}
};

// **********************************************
const int MAX = 1e5 + 5;
const int INF = 0x3f3f3f3f;

struct edge { int v, w; edge(int vv = 0, int ww = 0) : v(vv), w(ww) {} };

vector<edge> G[MAX];
int dis[MAX], vis[MAX];

void init(int n) {
	for (int i = 0; i <= n; i++) {
		G[i].clear();
		dis[i] = INF, vis[i] = 0;
	}
}

void addedge(int u, int v, int w) {
	G[u].emplace_back(v, w);
	G[v].emplace_back(u, w);
}

inline bool chkmin(int& x, int y) { return (x > y) ? (x = y, 1) : 0; }

void Dijk(int st, int ed) {
	priority_queue<pair<int, int>> Q;
	Q.emplace(0, st); dis[st] = 0;
	while (!Q.empty()) {
		int u = Q.top().second; Q.pop();
		if (u == ed) return;
		if (vis[u]) continue; else vis[u] = 1;
		for (edge& e : G[u]) {
			if (!vis[e.v] && chkmin(dis[e.v], dis[u] + e.w)) {
				Q.emplace(-dis[e.v], e.v);
			}
		}
	}
}