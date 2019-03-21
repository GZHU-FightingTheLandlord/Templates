struct biCut {
  vector<int> G[maxn];
  int N, tag, dfn[maxn], low[maxn];
  bool iscut[maxn];

  void init(int n) {
    N = n, tag = 0;
    for (int i = 1; i <= N; i++) {
      dfn[i] = low[i] = 0, iscut[i] = false;
      G[i].clear();
    }
  }
  void addedge(int u, int v) {
    G[u].push_back(v), G[v].push_back(u);
  }
  void dfs(int u, int f) {
    low[u] = dfn[u] = ++tag;
    int child = 0;
    for (auto& v : G[u]) {
      if (!dfn[v]) {
        ++child, dfs(v, u);
        low[u] = min(low[u], low[v]);
        if (low[v] >= dfn[u]) {
          iscut[u] = true;
        }
      } else if (dfn[v] < dfn[u] && v != f) {
        low[u] = min(low[u], dfn[v]);
      }
    }
    if (f < 0 && child == 1) {
      iscut[u] = false;
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!dfn[i]) dfs(i, -1);
    }
  }
};