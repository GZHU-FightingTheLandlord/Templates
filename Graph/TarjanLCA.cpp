const int N = 1e5 + 5;
const int M = 1e5 + 5;

vector<int> G[N];
vector<pair<int, int>> Q[N];
int lca[M], anc[N];
bool vis[N];

void init(int n) {
  iota(anc, anc + n, 0);
}

void addEdge(int u, int v);

void addQuery(int u, int v, int i) {
  Q[u].push_back({ v, i }), Q[v].push_back({ u, i });
}

int find(int x) {
  return x == anc[x] ? x : anc[x] = find(anc[x]);
}

void dfs(int u, int f) {
  for (auto& v : G[u]) {
    if (v != f) dfs(v, u), anc[v] = anc[u];
  }
  vis[u] = 1;
  for (auto& q : Q[u]) {
    int v = q.first, i = q.second;
    if (vis[v]) lca[i] = find(v);
  }
}