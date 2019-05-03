## Graph

### 2-sat

```cpp
struct twoSat {
  struct edge {
    int v, next;
    edge(int a = 0, int b = 0) : v(a), next(b) {}
  }G[maxm];
  int tot, head[maxn], mark[maxn], sz, stk[maxn];
  void init() {
    tot = 0;
    memset(mark, 0, sizeof mark);
    memset(head, -1, sizeof head);
  }
  // for every case u, (status[u] xor status[u ^ 1]) == true.
  // 
  // addcase: if status[u] == true then status[v] == true,
  // but if status[u] == false then status[v] can be true or false.
  // 
  void addcase(int u, int v) {
    G[tot] = edge(v, head[u]); head[u] = tot++;
  }
  int dfs(int u) {
    if (mark[u ^ 1]) return 0;
    if (mark[u]) return 1;
    stk[sz++] = u, mark[u] = 1;
    for (int i = head[u]; ~i; i = G[i].next) {
      if (!dfs(G[i].v)) return 0;
    }
    return 1;
  }
  int solve(int n) {
    for (int i = 0; i < n; i += 2) {
      if (!mark[i] && !mark[i ^ 1]) {
        sz = 0;
        if (!dfs(i)) {
          while (sz > 0) mark[stk[--sz]] = 0;
          if (!dfs(i ^ 1)) return 0;
        }
      }
    }
    return 1;
  }
}sat;
```

### 强连通

#### Tarjan

```cpp
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
```

#### kosaraju

```cpp
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
```

### 双连通

pending...

### 欧拉路

#### 无向

```cpp
// undirected, 0-base
template <int NV> class Hierholzer {
public:
  vector<int> path;
  multiset<int> G[NV];

  void addedge(int u, int v) {
    G[u].insert(v), G[v].insert(u);
  }

  void dfs(int cur) {
    while (!G[cur].empty()) {
      int tar = *G[cur].begin();
      G[cur].erase(G[cur].begin());
      G[tar].erase(G[tar].find(cur));
      dfs(tar);
    }
    path.push_back(cur);
  }

  bool get() {
    int src = -1, odd = 0, tot = 0;
    for (int i = 0; i < NV; i++) {
      tot += G[i].size();
      if (G[i].size() % 2 == 1) {
        odd++, src = (~src) ? src : i;
      }
    }
    if (odd != 0 && odd != 2) return false;
    dfs(odd ? src : 0);
    reverse(path.begin(), path.end());
    return (int)path.size() == tot / 2 + 1;
  }

  vector<int> get(int src) {
    dfs(src);
    reverse(path.begin(), path.end());
    return path;
  }
};
```

#### 有向

```cpp
// directed, 0-base.
template <int NV> class Hierholzer {
public:
  int deg[NV];
  vector<int> path;
  multiset<int> G[NV];

  void addedge(int u, int v) {
    G[u].insert(v), deg[u]++, deg[v]--;
  }

  void dfs(int cur) {
    while (!G[cur].empty()) {
      int tar = *G[cur].begin();
      G[cur].erase(G[cur].begin());
      dfs(tar);
    }
    path.push_back(cur);
  }

  bool get() {
    int src = -1, tot = 0, U = 0, D = 0, UZ = 0;
    for (int i = 0; i < NV; i++) {
      tot += G[i].size();
      if (deg[i] != 0) {
        U += (deg[i] == 1), D += (deg[i] == -1), UZ++;
        src = (~src) ? src : i;
      }
    }
    if (UZ != 0 && (UZ != 2 || U != 1 || D != 1)) return false;
    dfs(UZ ? src : 0);
    reverse(path.begin(), path.end());
    return (int)path.size() == tot + 1;
  }

  vector<int> get(int src) {
    dfs(src);
    reverse(path.begin(), path.end());
    return path;
  }
};
```

### 费用流

```cpp
struct edge {
  int v, cost, flow, cap, next;
  edge() {}
  edge(int V, int Cost, int Flow, int Cap, int nxt) : \
      v(V), cost(Cost), flow(Flow), cap(Cap), next(nxt) {}
} G[maxm << 1];
int tot, head[maxn], cost[maxn], inq[maxn], pre[maxn];

void init() {
  tot = 0;
  memset(head, -1, sizeof head);
}

void addedge(int u, int v, int cap, int cost) {
  G[tot] = edge(v, cost, 0, cap, head[u]); head[u] = tot++;
  G[tot] = edge(u, -cost, cap, cap, head[v]); head[v] = tot++;
}

bool spfa(int src, int dst) {
  memset(inq, 0, sizeof inq);
  memset(pre, -1, sizeof pre);
  memset(cost, 0x3f, sizeof cost);
  queue<int> Q; Q.push(src), cost[src] = 0;
  while (!Q.empty()) {
    int u = Q.front(); Q.pop(), inq[u] = 0;
    for (int i = head[u]; ~i; i = G[i].next) {
      edge &e = G[i];
      if (e.flow < e.cap && chkmin(cost[e.v], cost[u] + e.cost)) {
        pre[e.v] = i;
        if (!inq[e.v]) Q.push(e.v), inq[e.v] = 1;
      }
    }
  }
  return cost[dst] < 0x3f3f3f3f;
}

pair<int, int> mcmf(int src, int dst) {
  int totCost = 0, totFlow = 0;
  while (spfa(src, dst)) {
    int maxFlow = 0x3f3f3f3f;
    for (int u = dst; u != src; u = G[pre[u] ^ 1].v) {
      edge &e = G[pre[u]]; // , &r = G[pre[u] ^ 1];
      maxFlow = min(maxFlow, e.cap - e.flow);
    }
    totCost += maxFlow * cost[dst], totFlow += maxFlow;
    for (int u = dst; u != src; u = G[pre[u] ^ 1].v) {
      edge &e = G[pre[u]], &r = G[pre[u] ^ 1];
      e.flow += maxFlow, r.flow -= maxFlow;
    }
  }
  return { totFlow, totCost };
}
```

### 二分图匹配

```cpp
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
  int getans(int n) {
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
```

### 树剖(lca为例)

```cpp
int SZ[N], fa[N], son[N], top[N], dep[N];
int dfn, in[N], out[N];

void getsz(int u, int d, int f) {
  SZ[u] = 1, dep[u] = d, fa[u] = f;
  son[u] = 0;
  for (auto& v : G[u]) {
    if (v != f) {
      getsz(v, d + 1, u);
      SZ[u] += SZ[v];
      if (SZ[son[u]] < SZ[v]) son[u] = v;
    }
  }
}

void dfs(int u, int t) {
  in[u] = ++dfn, top[u] = t;
  if (son[u]) dfs(son[u], t);
  for (auto& v : G[u]) {
    if (v != fa[u] && v != son[u]) {
      dfs(v, v);
    }
  }
  out[u] = dfn;
}

int getlca(int u, int v) {
  for (; top[u] != top[v]; u = fa[top[u]]) {
    if (dep[top[u]] < dep[top[v]]) swap(u, v);
  }
  return dep[u] < dep[v] ? u : v;
}
```
