## DataStructure

### 区间增减树状数组

```cpp
struct Interval {
  int N, base[2][maxn];
  void setN(int n) { N = n; }
  void init() { memset(base, 0, sizeof base); }
  void add(int at, int v) {
    if (!at) return;
    for (int i = at; i <= N; i += i & -i) {
      base[0][i] += v, base[1][i] -= v * at;
    }
  }
  void add(int l, int r, int v) {
    add(l, v), add(r + 1, -v);
  }
  int getSum(int at) {
    int sum = 0, mul = at + 1;
    for (int i = at; i; i -= i & -i) {
      sum += mul * base[0][i] + base[1][i];
    }
    return sum;
  }
  int query(int l, int r) {
    return getSum(r) - getSum(l - 1);
  }
};
```

### 无旋Treap

```cpp
struct Treap {
#define ls(x) T[x].son[0]
#define rs(x) T[x].son[1]

  struct Node {
    int son[2], size, v, key, rev;
  } T[maxn];
  int tot, root;

  Treap() { tot = root = 0; }

  inline void init() { tot = root = 0; }
  
  inline void pushup(int i) {
    T[i].size = T[ls(i)].size + T[rs(i)].size + 1;
  }
  inline void pushdown(int i) {
    if (T[i].rev) {
      swap(ls(i), rs(i));
      T[ls(i)].rev ^= 1, T[rs(i)].rev ^= 1;
      T[i].rev = 0;
    }
  }
  void split(int rt, int &x, int &y, int v) {
    if (!rt) return (void)(x = y = 0);
    pushdown(rt);
    if (T[rt].v <= v) {
      x = rt, split(rs(rt), rs(x), y, v);
    } else {
      y = rt, split(ls(rt), x, ls(y), v);
    }
    pushup(rt);
  }
  void merge(int &rt, int x, int y) {
    if (!x || !y) {
      rt = x + y;
      return;
    }
    if (T[x].key < T[y].key) {
      pushdown(x), rt = x, merge(rs(rt), rs(x), y);
    } else {
      pushdown(y), rt = y, merge(ls(rt), x, ls(y));
    }
    pushup(rt);
  }
  inline void insert(int &rt, int v) {
    int x = 0, y = 0, z = ++tot;
    T[z].v = v, T[z].key = rand(), T[z].size = 1, T[z].rev = 0;
    split(rt, x, y, v), merge(x, x, z), merge(rt, x, y);
  }
  inline void erase(int &rt, int v) {
    int x = 0, y = 0, z = 0;
    split(rt, x, y, v), split(x, x, z, v - 1);
    merge(z, ls(z), rs(z)), merge(x, x, z), merge(rt, x, y);
  }
  inline int findkth(int rt, int k) {
    if (k == 0) return -inf;
    pushdown(rt);
    while (T[ls(rt)].size + 1 != k) {
      if (T[ls(rt)].size >= k) rt = ls(rt);
      else k -= (T[ls(rt)].size + 1), rt = rs(rt);
      pushdown(rt);
    }
    return T[rt].v;
  }
  inline int getrank(int &rt, int v) {
    int x = 0, y = 0, res;
    split(rt, x, y, v - 1), res = T[x].size + 1;
    return merge(rt, x, y), res;
  }
  inline int getpre(int &rt, int v) {
    int x = 0, y = 0, res;
    split(rt, x, y, v - 1), res = findkth(x, T[x].size);
    return merge(rt, x, y), res;
  }
  inline int getsuf(int &rt, int v) {
    int x = 0, y = 0, res;
    split(rt, x, y, v), res = findkth(y, 1);
    return merge(rt, x, y), res;
  }

  inline void insert(int v) { insert(root, v); }
  inline void erase(int v) { erase(root, v); }
  inline int findkth(int k) { return findkth(root, k); }
  inline int getrank(int v) { return getrank(root, v); }
  inline int getpre(int v) { return getpre(root, v); }
  inline int getsuf(int v) { return getsuf(root, v); }
} treap;
```

### ST表

```cpp
struct ST {
  vector<vector<int>> table;
  ST(vector<int> a = {}) {
    int n = a.size();
    table.resize(n, vector<int>(32 - __builtin_clz(n)));
    for (int i = 0; i < n; i++) {
      table[i][0] = a[i];
    }
    for (int j = 1; (1 << j) - 1 < n; j++) {
      for (int i = 0; i + (1 << j) - 1 < n; i++) {
        int x = table[i][j - 1], y = table[i + (1 << (j - 1))][j - 1];
        table[i][j] = min(x, y);
      }
    }
  }
  inline int getMin(int l, int r) {
    int k = 31 - __builtin_clz(r - l + 1);
    return min(table[l][k], table[r - (1 << k) + 1][k]);
  }
};
```

### 并查集(带权)

```cpp
template <int NV> class Dsu {
  int anc[NV], weight[NV];
  void init(int n = NV) {
    iota(anc.begin(), next(anc.begin(), n), 0);
    fill(anc.begin(), next(anc.begin(), n), 0);
  }
  int find(int x) {
    if (x == anc[x]) return x;
    int fa = anc[x];
    anc[x] = find(anc[x]);
    weight[x] += weight[fa];
    return anc[x];
  }
  bool unite(int u, int v, int w = 0) {
    int a = find(u), b = find(v);
    if (a == b) return false;
    anc[b] = a;
    weight[b] = weight[u] + w - weight[v];
    return true;
  }
};
```