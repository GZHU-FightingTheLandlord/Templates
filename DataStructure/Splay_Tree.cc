// Example: https://ac.nowcoder.com/acm/contest/141/C

#include <stdio.h>
#include <algorithm>
using std::swap;

struct Splay_Tree {
  static const int maxn = 2e5 + 5;
  static const int INF = 0x3f3f3f3f;

  struct Node {
    int val, rev;
    int son[2], size;
    void init(int _v) {
      val = _v, size = 1;
      rev = son[0] = son[1] = 0;
    }
  }T[maxn];
  int fa[maxn], root;

  void init(int n) {
    T[0].init(-INF); T[1].init(-INF);  T[n + 2].init(-INF);
    // set as [1, n]
    for (int i = 2; i <= n + 1; i++) {
      T[i].init(i - 1);
    }
    root = build(1, n + 2);
    fa[0] = 0; T[0].son[1] = root; T[0].size = 0;
  }

  int build(int l, int r) {
    if (l > r) return 0;
    if (l == r) return l;
    int mid = (l + r) >> 1, sl, sr;
    T[mid].son[0] = sl = build(l, mid - 1);
    T[mid].son[1] = sr = build(mid + 1, r);
    fa[sl] = fa[sr] = mid;
    pushup(mid);
    return mid;
  }

  void pushup(int x) {
    T[x].size = 1;
    if (T[x].son[0]) {
      T[x].size += T[T[x].son[0]].size;
    }
    if (T[x].son[1]) {
      T[x].size += T[T[x].son[1]].size;
    }
  }

  void pushdown(int x) {
    if (x == 0) return;
    if (T[x].rev) {
      if (T[x].son[0]) T[T[x].son[0]].rev ^= 1;
      if (T[x].son[1]) T[T[x].son[1]].rev ^= 1;
      swap(T[x].son[0], T[x].son[1]);
      T[x].rev = 0;
    }
  }

  // zig or zag
  void rotate(int x, int kind) {
    int y = fa[x], z = fa[y];
    T[y].son[!kind] = T[x].son[kind]; fa[T[x].son[kind]] = y;
    T[x].son[kind] = y; fa[y] = x;
    T[z].son[T[z].son[1] == y] = x; fa[x] = z;
    pushup(y);
  }

  // let x become tar's child
  // when tar == 0, x will be the root
  void Splay(int x, int tar) {
    if (x == tar) return;
    while (fa[x] != tar) {
      int y = fa[x], z = fa[y];
      pushdown(z); pushdown(y); pushdown(x);
      int rx = T[y].son[0] == x, ry = T[z].son[0] == y;
      if (z == tar) rotate(x, rx);
      else {
        if (rx == ry) rotate(y, ry);
        else rotate(x, rx);
        rotate(x, ry);
      }
    }
    pushup(x);
    if (tar == 0) root = x;
  }

  // get the kth Node
  int findkth(int pos) {
    int u = root;
    pushdown(u);
    while (T[T[u].son[0]].size != pos) {
      if (pos < T[T[u].son[0]].size) u = T[u].son[0];
      else {
        pos -= T[T[u].son[0]].size + 1;
        u = T[u].son[1];
      }
      pushdown(u);
    }
    return u;
  }

  // pre-order travel
  void print(int rt) {
    pushdown(rt);
    if (T[rt].son[0]) print(T[rt].son[0]);
    if (T[rt].val != -INF) printf("%d ", T[rt].val);
    if (T[rt].son[1]) print(T[rt].son[1]);
  }

  // reverse [l, r]
  void reverse(int l, int r) {
    int u = findkth(l - 1), v = findkth(r + 1);
    Splay(u, 0); Splay(v, u);
    T[T[v].son[0]].rev ^= 1;
  }

  // do sth.
}ST;
