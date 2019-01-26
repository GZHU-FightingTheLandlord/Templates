const int maxn = 1e5 + 5;
const int maxm = 2e5 + 5;

struct edge {
  int v, next;
  edge(int a = 0, int b = 0) : v(a), next(b) {}
};

struct Scc {
  edge G[maxm];
  int tot, head[maxn];

  int index, count;
  int stack[maxn], top;
  int dfn[maxn], low[maxn], ins[maxn], sccno[maxn];

  void init() {
    tot = 0;
    memset(head, -1, sizeof head);
    memset(sccno, -1, sizeof sccno);
  }

  void addedge(int u, int v) {
    G[tot] = edge(v, head[u]); head[u] = tot++;
    G[tot] = edge(u, head[v]); head[v] = tot++;
  }

  void dfs(int u) {
    dfn[u] = low[u] = ++index;
    stack[top++] = u, ins[u] = 1;

    int v;
    for (int i = head[u]; ~i; i = G[i].next) {
      v = G[i].v;
      if (!dfn[v]) {
        dfs(v);
        low[u] = min(low[u], low[v]);
      } else if (ins[v]) {
        low[u] = min(low[u], dfn[v]);
      }
    }
    if (dfn[u] == low[u]) {
      ++count;
      do {
        v = stack[--top];
        ins[v] = 0, sccno[v] = count;
      } while (u != v);
    }
  }

  void solve(int n) {
    index = count = top = 0;
    memset(dfn, 0, sizeof dfn);
    memset(low, 0, sizeof low);
    
    for (int i = 1; i <= n; i++) {
      if (!dfn[i]) {
        dfs(i);
      }
    }
  }

  int operator[] (int x) const { return sccno[x]; }
};