#include <algorithm>
#include <vector>
#include <queue>
using namespace std;

struct MCMF {
    struct edge {
        int u, v, c, w;
        edge(int uu = 0, int vv = 0, int cc = 0, int ww = 0) {
            u = uu, v = vv, c = cc, w = ww;
        }
    };

    int N, cnt, st, ed, SumCost, SumFlow;
    vector<edge> E;
    vector<vector<int>> G;
    vector<int> cost, pre, inq;

    MCMF(int n = 0, int ecnt = 0) : N(n), E(n + 5), G(n + 5), cost(n + 5),\
    inq(n + 5), pre(n + 5) { init(); }

    void init(int n = -1)
    {
        if (n != -1) N = n;

        cnt = 0;
        for (int i = 0; i <= N; i++) {
            G[i] = vector<int>();
        }
    }

    void addedge(int u, int v, int c, int w)
    {
        E[cnt] = edge(u, v, c, w); G[u].push_back(cnt); cnt++;
        E[cnt] = edge(v, u, -c, 0); G[v].push_back(cnt); cnt++;
    }

    bool spfa(int st, int ed)
    {
        fill(cost.begin(), cost.end(), 0x3f3f3f3f);
        fill(inq.begin(), inq.end(), 0);
        fill(pre.begin(), pre.end(), -1);
        queue<int> Q;
        Q.push(st); inq[st] = 1; cost[st] = 0;
        while (!Q.empty()) {
            int u = Q.front(); Q.pop(); inq[u] = 0;

            for (int i = 0; i < G[u].size(); i++) {
                int v = E[G[u][i]].v, c = E[G[u][i]].c;
                if (E[G[u][i]].w > 0 && cost[v] > cost[u] + c) {
                    cost[v] = cost[u] + c;
                    pre[v] = G[u][i];
                    if (!inq[v]) {
                        Q.push(v);
                        inq[v] = 1;
                    }
                }
            }
        }
        return pre[ed] != -1;
    }

    pair<int, int> augment(int st, int ed)
    {
        int ret = 0x3f3f3f3f;
        for (int now = ed, edg; now != st; now = E[pre[now] ^ 1].v) {
            edg = pre[now];
            ret = min(ret, E[edg].v);
        }
        for (int now = ed, edg; now != st; now = E[pre[now] ^ 1].v) {
            edg = pre[now];
            E[edg].w -= ret;
            E[edg ^ 1].w += ret;
        }
        return make_pair(cost[ed], ret);
    }

    void solve(int st, int ed)
    {
        SumCost = 0, SumFlow = 0;
        while (spfa(st, ed)) {
            pair<int, int> tmp = augment(st, ed);
            SumCost += tmp.first, SumFlow += tmp.second;
        }
    }
};


// https://vjudge.net/problem/HDU-6611
// dijkstra MCMF
const int INF = 0x3f3f3f3f;
namespace MCMF {
  struct edge {
    int to, cap, cost, rev;
    edge() {}
    edge(int tt, int ca, int co, int re) {
      to = tt, cap = ca, cost = co, rev = re;
    }
  };
  const int N = 5000;
  int V, H[N + 5], dis[N + 5], prev[N + 5], pree[N + 5];
  vector<edge> G[N + 5];
  void init(int n) {
    V = n;
    for(int i = 0; i <= V; i++) {
      G[i].clear();
    }
  }
  void addedge(int from, int to, int cap, int cost) {
    G[from].push_back(edge(to, cap, cost, G[to].size()));
    G[to].push_back(edge(from, 0, -cost, G[from].size() - 1));
  }
  // what is f and flow = =
  // int solve(int s, int t, int f, int &flow) { 
  int solve(int s, int t, int f, int flow) {
    int res = 0;
    memset(H, 0, sizeof(*H) * (V + 1));
    while(f) {
      priority_queue<pii, vector<pii>, greater<pii>> q;
      memset(dis, 0x3f, sizeof(*dis) * (V + 1));
      dis[s] = 0;
      q.push({0, s});
      while(!q.empty()) {
        pii now = q.top();
        q.pop();
        int v = now.second;
        if(dis[v] < now.first) {
          continue;
        }
        for(int i = 0; i < int(G[v].size()); i++) {
          edge &e = G[v][i];
          if(e.cap > 0 && dis[e.to] > dis[v] + e.cost + H[v] - H[e.to]) {
            dis[e.to] = dis[v] + e.cost + H[v] - H[e.to];
            prev[e.to] = v;
            pree[e.to] = i;
            q.push({dis[e.to], e.to});
          }
        }
      }
      if(dis[t] == INF) break;
      for(int i = 0; i <= V; i++) {
        H[i] += dis[i];
      }
      int d = f;
      for(int v = t; v != s; v = prev[v]) {
        d = min(d, G[prev[v]][pree[v]].cap);
      }
      f -= d; flow += d; res += d * H[t];
      for(int v = t; v != s; v = prev[v]) {
        edge &e = G[prev[v]][pree[v]];
        e.cap -= d;
        G[v][e.rev].cap += d;
      }
    }
    return res;
  }
}
