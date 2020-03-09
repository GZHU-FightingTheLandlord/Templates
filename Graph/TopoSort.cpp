template <int NV> struct bfs {
  std::vector<int> G[NV], ord;
  int deg[NV];

  bfs(int n = 0) {
    init(n);
  }
  void init(int n) {
    ord.clear();
    for (int i = 0; i < n; i++) {
      G[i].clear(), deg[i] = 0;
    }
  }
  void addedge(int u, int v) {
    G[u].push_back(v), ++deg[v];
  }
  void gao(int n) {
    std::queue<int> Q;
    for (int i = 0; i < n; i++) {
      if (!deg[i]) {
        Q.push(i);
      }
    }
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      ord.push_back(u);
      for (auto v : G[u]) {
        if (--deg[v] == 0) {
          Q.push(v);
        }
      }
    }
  }
};
