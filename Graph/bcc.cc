// vertex-bcc
struct bcc {
  struct edge { int u, v; };
  vector<int> G[maxn], cont[maxn];
  int N, tag, tot, dfn[maxn], bccno[maxn];
  bool iscut[maxn];
  stack<edge> S;

  void init(int n) {
    N = n, tag = tot = 0;
    for (int i = 1; i <= N; i++) {
      G[i].clear();
      dfn[i] = bccno[i] = 0;
      iscut[i] = false;
    }
  }
  void addedge(int u, int v) {
    G[u].push_back(v), G[v].push_back(u);
  }
  int dfs(int u, int f) {
    int lowu = dfn[u] = ++tag;
    int child = 0;
    for (auto& v : G[u]) {
      if (!dfn[v]) {
        ++child, S.push({ u, v });
        int lowv = dfs(v, u);
        lowu = min(lowu, lowv);
        if (lowv >= dfn[u]) {
          iscut[u] = true;
          cont[++tot].clear();
          while (true) {
            edge e = S.top(); S.pop();
            if (bccno[e.u] != tot) {
              cont[tot].push_back(e.u);
              bccno[e.u] = tot;
            }
            if (bccno[e.v] != tot) {
              cont[tot].push_back(e.v);
              bccno[e.v] = tot;
            }
            if (e.u == u && e.v == v) {
              break;
            }
          }
        }
      } else if (dfn[v] < dfn[u] && v != f) {
        S.push({ u, v });
        lowu = min(lowu, dfn[v]);
      }
    }
    if (f < 0 && child == 1) {
      iscut[u] = false;
    }
    return lowu;
  }
};