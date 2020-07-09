// 0-base, undirected
class hierholzer {
public:
  vector<int> path;
  vector<multiset<int>> G;

  hierholzer(int n) {
    init(n);
  }

  void init(int n) {
    path.clear();
    G.assign(n, multiset<int>());
  }
  
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
    for (int i = 0; i < (int) G.size(); i++) {
      tot += G[i].size();
      if (G[i].size() % 2 == 1) {
        odd++, src = (~src) ? src : i;
      }
    }
    if (odd != 0 && odd != 2) return false;
    dfs(odd ? src : 0);
    reverse(path.begin(), path.end());
    return (int) path.size() == tot / 2 + 1;
  }

  vector<int> get(int src) {
    dfs(src);
    reverse(path.begin(), path.end());
    return path;
  }
};

// directed, 0-base.
template <int NV> class HierholzerD {
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