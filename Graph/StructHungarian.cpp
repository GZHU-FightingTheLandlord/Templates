struct maxMatch {
  int link[maxn], vis[maxn];
  bool find(int u) {
    for (int i = head[u]; i != -1; i = G[i].next) {
      int v = G[i].v;
      if (!vis[v]) {
        vis[v] = 1;
        if (link[v] == -1 || find(link[v])) {
          link[v] = u;
          // link[u] = v;
          return true;
        }
      }
    }
    return false;
  }
  int getAns(int n) {
    int ans = 0;
    memset(link, -1, sizeof link);
    for (int i = 1; i <= n; i++) {
      if (link[i] == -1) {
        memset(vis, 0, sizeof vis);
        if (find(i)) ++ans;
      }
    }
    return ans;
  }
};