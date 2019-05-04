struct edge {
  int v, next;
} G[M];
int tot, h[N], ord, dfn[N], low[N];
bool iscut[N], isbridge[M];

void init() {
  tot = ord = 0;
  memset(h, -1, sizeof h);
  memset(dfn, 0, sizeof dfn);
  memset(low, 0, sizeof low);
  memset(iscut, false, sizeof iscut);
  memset(isbridge, false, sizeof false);
}

void addedge(int u, int v) {
  G[tot] = { v, h[u] }, h[u] = tot++;
  G[tot] = { u, h[v] }, h[v] = tot++;
}

void dfs(int u, int f) {
  low[u] = dfn[u] = ++ord;
  int child = 0;
  for (int i = h[u]; ~i; i = G[i].next) {
    edge &e = G[i];
    if (!dfn[e.v]) {
      ++child, dfs(e.v, u);
      low[u] = min(low[u], low[e.v]);
      if (low[e.v] >= dfn[u]) {
        iscut[u] = true;
      }
      if (low[e.v] > dfn[u]) {
        isbridge[i] = isbridge[i ^ 1] = true;
      }
    } else if (dfn[e.v] < dfn[u] && e.v != f) {
      low[u] = min(low[u], dfn[e.v]);
    }
  }
  if (f == -1 && child == 1) {
    iscut[u] = false;
  }
}

void solve(int n) {
  for (int i = 1; i <= n; i++) {
    if (!dfn[i]) dfs(i, -1);
  }
}