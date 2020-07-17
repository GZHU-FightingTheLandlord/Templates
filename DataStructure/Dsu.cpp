struct Dsu {
  vector<int> anc, size;
  Dsu(int n = 0) : anc(n), size(n, 0) {
    iota(anc.begin(), anc.end(), 0);
  }
  int operator[] (int x) {
    return x == anc[x] ? x : anc[x] = operator[](anc[x]);
  }
  bool operator() (int u, int v) {
    int a = operator[](u), b = operator[](v);
    if (a == b) return false;
    if (size[a] < size[b]) anc[a] = b;
    else anc[b] = a, size[a] += (size[a] == size[b]);
    return true;
  }
};
