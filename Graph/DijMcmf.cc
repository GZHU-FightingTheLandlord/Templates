#include <bits/stdc++.h>
using namespace std;

const int INF = 0x3f3f3f3f;
struct mcmf {
    static const int maxn = 205;
    static const int maxm = 2e5 + 5;
    struct edge {
        int v, w, c, next;
        edge(int a = 0, int b = 0, int d = 0, int e = 0) : v(a), w(b), c(d), next(e) {}
    }G[maxm];

    int tot, n, head[maxn], cost[maxn], h[maxn], pre[maxn];
    bool vis[maxn];
    
    mcmf() {}
    
    inline void init(int N) { n = N, tot = 0, memset(head, -1, sizeof head); }

    inline void addedge(int u, int v, int w, int c) {
        G[tot] = edge(v, w, c, head[u]); head[u] = tot++;
        G[tot] = edge(u, 0, -c, head[v]); head[v] = tot++;
    }

    inline void spfa(int s, int t) {
        queue<int> Q;
        memset(h, 0x3f, sizeof h);
        memset(vis, false, sizeof vis);
        Q.push(s), h[s] = 0, vis[s] = true;
        while (!Q.empty()) {
            int u = Q.front(); Q.pop(), vis[u] = false;
            for (int i = head[u]; i != -1; i = G[i].next) {
                edge &e = G[i];
                if (e.w > 0 && h[e.v] > h[u] + e.c) {
                    h[e.v] = h[u] + e.c;
                    if (!vis[e.v]) Q.push(e.v), vis[e.v] = true;
                }
            }
        }
    }

    inline bool Dijstra(int s, int t) {
        memset(pre, -1, sizeof pre);
        memset(vis, false, sizeof vis);
        memset(cost, 0x3f, sizeof cost);
        priority_queue<pair<int, int>> Q;
        Q.emplace(0, s), cost[s] = 0;
        while (!Q.empty()) {
            int u = Q.top().second; Q.pop();
            if (u == t) return true;
            if (vis[u]) continue;
            vis[u] = true;
            for (int i = head[u]; i != -1; i = G[i].next) {
                int v = G[i].v, c = G[i].c + h[u] - h[v];
                if (G[i].w > 0 && !vis[v] && cost[v] > cost[u] + c) {
                    cost[v] = cost[u] + c;
                    pre[v] = i;
                    Q.emplace(-cost[v], v);
                }
            }
        }
        return false;
    }

    int solve(int s, int t) {
        spfa(s, t);
        int res = 0, cst = 0;
        while (Dijstra(s, t)) {
            int minflow = INF;
            for (int i = t; i != s; i = G[pre[i] ^ 1].v) {
                minflow = min(minflow, G[pre[i]].w);
            }
            for (int i = t; i != s; i = G[pre[i] ^ 1].v) {
                int j = pre[i];
                G[j].w -= minflow;
                G[j ^ 1].w += minflow;
            }
            res += minflow;
            cst += minflow * (h[t] + cost[t]);
            for (int i = 1; i <= n; i++) {
                h[i] = min(h[i] + cost[i], INF);
            }
        }
        return cst;
    }
}M;