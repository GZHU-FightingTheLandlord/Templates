struct biBridge {
  vector<pair<int, int>> G[maxn];
  int N, tag, dfn[maxn], low[maxn];
  bool isbridge[maxm];

  void init(int n) {
    N = n, tag = 0;
    for (int i = 1; i <= n; i++) {
      dfn[i] = low[i] = 0, G[i].clear();
    }
    memset(isbridge, 0, sizeof isbridge);
  }
  void addedge(int u, int v, int i) {
    G[u].push_back({ v, i });
    G[v].push_back({ u, i });
  }
  void dfs(int u, int f) {
    dfn[u] = low[u] = ++tag;
    int child = 0;
    for (auto& e : G[u]) {
      int v = e.first, index = e.second;
      if (!dfn[v]) {
        ++child, dfs(v, u);
        low[u] = min(low[u], low[v]);
        if (low[v] > dfn[u]) {
          isbridge[index] = true;
        }
      } else if (dfn[v] < dfn[u] && v != f) {
        low[u] = min(low[u], dfn[v]);
      }
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!dfn[i]) dfs(i, -1);
    }
  }
} bridge;