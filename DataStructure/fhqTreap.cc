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