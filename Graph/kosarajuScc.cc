struct kosaraju {
  int N, tot, scc[maxn], vis[maxn];
  vector<int> G[maxn], R[maxn], acc;
  void init(int n) {
    N = n;
    tot = 0, acc.clear();
    for (int i = 1; i <= N; i++) {
      G[i].clear(), R[i].clear();
      vis[i] = 0, scc[i] = 0;
    }
  }
  void DFS1(int u) {
    vis[u] = 1;
    for (auto& v : G[u]) {
      if (!vis[v]) DFS1(v);
    }
    acc.push_back(u);
  }
  void DFS2(int u, int p) {
    scc[u] = p;
    for (auto& v : R[u]) {
      if (!scc[v]) DFS2(v, p);
    }
  }
  void solve() {
    for (int i = 1; i <= N; i++) {
      if (!vis[i]) DFS1(i);
    }
    reverse(acc.begin(), acc.end());
    for (auto& u : acc) {
      if (!scc[u]) DFS2(u, ++tot);
    }
  }
};