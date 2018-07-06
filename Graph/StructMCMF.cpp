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