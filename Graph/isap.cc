const int N = 1e2 + 5;
const int M = 2e4 + 5;
const int inf = 0x3f3f3f3f;

struct edge {
  int v, flow, cap, next;
} G[M];
int tot, n, src, dst, h[N], cur[N], gap[N], dep[N];

void init() {
  tot = 0;
  memset(h, -1, sizeof h);
}

void addedge(int u, int v, int w) {
  G[tot] = { v, 0, w, h[u] }, h[u] = tot++;
  G[tot] = { u, w, w, h[v] }, h[v] = tot++;
}

void bfs() {
  memset(gap, 0, sizeof gap);
  memset(dep, -1, sizeof dep);

  queue<int> Q; Q.push(dst);
  dep[dst] = 0, gap[0] = 1;
  while (!Q.empty()) {
    int u = Q.front(); Q.pop();
    for (int i = h[u]; ~i; i = G[i].next) {
      int v = G[i].v;
      if (~dep[v]) continue;
      Q.push(v), dep[v] = dep[u] + 1, gap[dep[v]]++;
    }
  }
}

int dfs(int u, int flow) {
  if (u == dst) return flow;
  int used = 0;
  for (int &i = cur[u]; ~i; i = G[i].next) {
    edge &e = G[i];
    if (e.flow < e.cap && dep[e.v] + 1 == dep[u]) {
      int tmp = dfs(e.v, min(e.cap - e.flow, flow - used));
      if (tmp == 0) continue;
      e.flow += tmp, G[i ^ 1].flow -= tmp, used += tmp;
      if (used == flow) return used;
    }
  }
  --gap[dep[u]];
  if (!gap[dep[u]]) dep[src] = n + 1;
  ++gap[++dep[u]];
  return used;
}

int isap() {
  bfs();
  int res = 0;
  while (dep[src] < n) {
    memcpy(cur, h, sizeof h);
    res += dfs(src, inf);
  }
  return res;
}