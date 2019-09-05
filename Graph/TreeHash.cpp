// 树哈希 预处理O(nln(n))（线性筛） 对每一个图 O(n+m)
// 标号从1开始 先使用init()初始化
// 然后对于每一个图先clear(n)再addedge(u,v)再solve()得出结果
namespace TH {
  typedef pair<int, int> pii;
  const int N = 1e6 + 5, M = 1e6 + 5, PN = 1.6e7 + 10;
  struct edge {
    int v, next;
  } G[M];
  int tot, h[N], siz[N], Siz, n;
  void clear(int _n) {
    tot = 0;
    n = _n;
    memset(h, 0xff, sizeof(*h) * (n + 1));
  }
  vector<int> pri;
  bool isnp[PN];
  void init() {
    isnp[0] = isnp[1] = true;
    for(int i = 2; i < PN; i++) {
      if(!isnp[i]) pri.push_back(i);
      for(int j = 0; j < int(pri.size()) && i * pri[j] < PN; j++) {
        isnp[i * pri[j]] = true;
        if(i % pri[j] == 0) break;
      }
    }
  }
  void addedge(int u, int v) {
    G[tot] = {v, h[u]}, h[u] = tot++;
    G[tot] = {u, h[v]}, h[v] = tot++;
  }
  pii root;
  void dfs(int u, int fa) {
    siz[u] = 1;
    int res = 0;
    for(int i = h[u]; ~i; i = G[i].next) {
      if(G[i].v == fa) continue;
      dfs(G[i].v, u);
      siz[u] += siz[G[i].v];
      res = max(res, siz[G[i].v]);
    }
    res = max(res, n - siz[u]);
    if(res < Siz) {
      root = {u, -1};
      Siz = res;
    } else if(res == Siz) {
      root.second = u;
    }
  }
  void getRoot() {
    Siz = 0x3f3f3f3f;
    root = {-1, -1};
    dfs(1, -1);
  }
  size_t Dfs(int u, int fa) {
    siz[u] = 1;
    size_t ret = 1;
    for(int i = h[u]; ~i; i = G[i].next) {
      if(G[i].v == fa) continue;
      size_t cur = Dfs(G[i].v, u);
      siz[u] += siz[G[i].v];
      cur *= size_t(pri[siz[G[i].v]]);
      ret += cur;
    }
    return ret;
  }
  pair<size_t, size_t> solve() {
    getRoot();
    pair<size_t, size_t> ret = {Dfs(root.first, -1), -1};
    if(root.second != -1) ret.second = Dfs(root.second, -1);
    if(ret.first > ret.second) swap(ret.first, ret.second);
    return ret;
  }
}
