struct Scc {
  vector<int> G[maxn];
  int N, tag, tot, dfn[maxn], low[maxn], sccno[maxn];
  stack<int> S;

  void init(int n) {
    N = n, tag = tot = 0;
    for (int i = 1; i <= n; i++) {
      dfn[i] = low[i] = sccno[i] = 0;
      G[i].clear();
    }
  }
  void addedge(int u, int v) {
    G[u].push_back(v);
  }
  void dfs(int u) {
    dfn[u] = low[u] = ++tag;
    S.push(u);
    for (auto& v : G[u]) {
      if (!dfn[v]) {
        dfs(v);
        low[u] = min(low[u], low[v]);
      } else if (!sccno[v]) {
        low[u] = min(low[u], dfn[v]);
      }
    }
    if (low[u] == dfn[u]) {
      ++tot;
      while (true) {
        int x = S.top(); S.pop();
        sccno[x] = tot;
        if (x == u) break;
      }
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!dfn[i]) dfs(i);
    }
  }
} scc;