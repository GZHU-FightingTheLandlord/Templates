# title

## link-cut tree

```cpp
namespace lct {
  struct Node {
    int size, fa, son[2];
    bool rev;
  } tree[N];

  bool isRoot(int x) {
    return tree[tree[x].fa].son[0] != x && tree[tree[x].fa].son[1] != x;
  }

  int which(int x) {
    return tree[tree[x].fa].son[1] == x;
  }

  void apply(int x) {
    swap(tree[x].son[0], tree[x].son[1]);
    tree[x].rev ^= 1;
  }

  void pushDown(int x) {
    if (tree[x].rev) {
      apply(tree[x].son[0]), apply(tree[x].son[1]);
      tree[x].rev = 0;
    }
  }

  void pushUp(int x) {
    tree[x].size = 1 + tree[tree[x].son[0]].size + tree[tree[x].son[1]].size;
  }

  void rotate(int x) {
    int y = tree[x].fa, z = tree[y].fa;
    int id = which(x), p = tree[x].son[1 - id];
    if (!isRoot(y)) tree[z].son[which(y)] = x;
    tree[x].fa = z, tree[y].son[id] = p, tree[p].fa = y;
    tree[x].son[1 - id] = y, tree[y].fa = x;
    pushUp(y), pushUp(x);
  }

  void dfs(int root) {
    if (!isRoot(root)) dfs(tree[root].fa);
    pushDown(root);
  }

  void Splay(int x) {
    dfs(x);
    while (!isRoot(x)) {
      int y = tree[x].fa;
      if (!isRoot(y)) {
        rotate((which(x) == which(y)) ? y : x);
      }
      rotate(x);
    }
  }

  // access u -> root
  int access(int u) {
    int tmp = 0;
    while (u != 0) {
      Splay(u), tree[u].son[1] = tmp, pushUp(u);
      tmp = u, u = tree[u].fa;
    }
    return tmp;
  }

  // find the root of this tree.
  int findRoot(int root) {
    access(root), Splay(root);
    int tmp = root;
    while (tree[tmp].son[0]) {
      pushDown(tmp);
      tmp = tree[tmp].son[0];
    }
    Splay(tmp);
    return tmp;
  }

  // make the node become root.
  void makeRoot(int root) {
    access(root), Splay(root), apply(root);
  }

  // link edge u -> v
  void link(int u, int v) {
    makeRoot(u); tree[u].fa = v;
  }

  // cut edge u -> v
  void cut(int u, int v) {
    makeRoot(u), access(v), Splay(v);
    tree[v].son[0] = tree[u].fa = 0, pushUp(v);
  }
}
```

## Steiner Tree

$ O(n3^k) $

```cpp
// maximum node num, key node num, infinity.
const int N = 1010, K = 5, inf = 0x3f3f3f3f;

struct edge { int v, w; };

vector<edge> G[N];
int dp[1 << K][N], a[N];
bool inq[N];
queue<int> Q;

void init() {
  memset(dp, 0x3f, sizeof dp);
  for (int i = 0; i < N; i++) {
    G[i].clear();
  }
}

void addEdge(int u, int v, int w);

void bellmanFord(int *dp) {
  while (!Q.empty()) {
    int u = Q.front(); Q.pop(), inq[u] = 0;
    for (auto& e : G[u]) {
      if (dp[e.v] > dp[u] + e.w) {
        dp[e.v] = dp[u] + e.w;
        if (!inq[e.v]) {
          Q.push(e.v), inq[e.v] = 1;
        }
      }
    }
  }
}

// k -> key node num, n -> node num.
void solve(int k, int n) {
  for (int S = 1; S < (1 << k); S++) {
    for (int i = 0; i < n; i++) {
      for (int s = (S - 1) & S; s; s = (s - 1) & S) {
        dp[S][i] = min(dp[S][i], dp[s][i] + dp[s ^ S][i]);
      }
      if (dp[S][i] < inf) Q.push(i), inq[i] = 1;
    }
    bellmanFord(dp[S]);
  }
}
```

## ZhuLiu's Algo

$O(nm)$(?)

```cpp
const int N = 1010, M = 1010101, inf = 0x3f3f3f3f;

// directed edges.
// r -> root, x_i -> u, y_i -> v, z_i -> w
int r, dfn, x[M], y[M], z[M], last[N], weight[N], id[N];
bool vis[N], instk[N];

void dfs(int u) {
  if (u == r) return;
  instk[u] = vis[u] = 1;
  if (!vis[last[u]]) {
    dfs(last[u]);
  } else if (instk[last[u]]) {
    id[u] = ++dfn;
    for (int v = last[u]; v != u; v = last[v]) {
      id[v] = dfn;
    }
  }
  instk[u] = 0;
  if (!id[u]) id[u] = ++dfn;
}

// n -> nodes num, m -> edges num.
int solve(int n, int m) {
  int ans = 0;
  while (true) {
    for (int i = 1; i <= n; i++) {
      weight[i] = inf;
    }
    for (int i = 1; i <= m; i++) {
      if (z[i] < weight[y[i]]) {
        last[y[i]] = x[i], weight[y[i]] = z[i];
      }
    }
    memset(id, 0, sizeof id);
    memset(vis, 0, sizeof vis);
    id[r] = dfn = 1;
    for (int i = 1; i <= n; i++) {
      if (i == r) continue;
      if (weight[i] == inf) return -1;
      ans += weight[i];
      if (!vis[i]) dfs(i);
    }
    if (dfn == n) return ans;
    int cnt = 0;
    for (int i = 1; i <= m; i++) {
      if (id[x[i]] != id[y[i]]) {
        z[++cnt] = z[i] - weight[y[i]];
        x[cnt] = id[x[i]], y[cnt] = id[y[i]];
      }
    }
    m = cnt, n = dfn, r = id[r];
  }
  return -1;
}
```

## CartesianTree

```cpp
// left son, right son, father.
int l[N], r[N], fa[N];

void build(int n, int arr[]) {
  static int top, stack[N];
  top = 0, stack[top++] = 1;
  for (int i = 2; i <= n; i++) {
    while (top && arr[stack[top - 1]] > arr[i]) --top;
    if (top > 0) {
      int x = i, y = stack[top - 1];
      l[x] = r[y], fa[l[x]] = x, fa[x] = y, r[y] = x;
      stack[top++] = x;
    } else {
      fa[stack[0]] = i, l[i] = stack[0];
      stack[top++] = i;
    }
  }  
}
```
