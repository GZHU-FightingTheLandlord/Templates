#include <algorithm>
using namespace std;

const int maxn = 2e5 + 5;

struct Node {
  int val, max, size;
  Node *fa, *ch[2];
  Node(int v = 0, Node *pre = NULL) : val(v), fa(pre) {
    ch[0] = ch[1] = NULL;
  }
  void pushup() {
    max = val, size = 1;
    if (ch[0] != NULL) size += ch[0]->size, max = std::max(max, ch[0]->max);
    if (ch[1] != NULL) size += ch[1]->size, max = std::max(max, ch[1]->max);
  }
  void pushdown() {
    // null
  }
};

struct SplayTree {
  int count;
  Node *root, box[maxn];

  void init() {
    root = NULL, count = 0;
  }
  // arr[0]与arr[n+1]需用作哨兵
  // build(ST.root, NULL, 0, n + 1, arr)
  void build(Node *&u, Node *fa, int l, int r, int *arr) {
    if (l > r) return (void)(u = NULL);
    int mid = (l + r) / 2;
    *(u = &box[++count]) = Node(arr[mid], fa);
    build(u->ch[0], u, l, mid - 1, arr);
    build(u->ch[1], u, mid + 1, r, arr);
    u->pushup();
  }
  // zig when c = 0, zag when c = 1
  void rotate(Node *x, int c) {
    Node *y = x->fa, *z = y->fa;
    y->ch[!c] = x->ch[c];
    if (x->ch[c]) x->ch[c]->fa = y;
    x->fa = z;
    if (z) z->ch[z->ch[1] == y] = x;
    y->fa = x, y->pushup(), x->ch[c] = y;
  }
  // 把x提到target下, 使x->fa == target
  // target = NULL时, 把x提到root
  void Splay(Node *x, Node *target) {
    if (x == target) return;
    while (x->fa != target) {
      Node *y = x->fa, *z = y->fa;
      if (z == target) rotate(x, y->ch[0] == x);
      else {
        if (z->ch[0] == y) {
          if (y->ch[0] == x) rotate(y, 1);
          else rotate(x, 0);
          rotate(x, 1);
        } else {
          if (y->ch[1] == x) rotate(y, 0);
          else rotate(x, 1);
          rotate(x, 0);
        }
      }
    }
    x->pushup();
    if (target == NULL) root = x;
  }
  // 把kth node提到target下
  // target == NULL时, 把kth node提到root
  void select(int k, Node *target) {
    Node *t = root;
    for (int sz = 0; true; ) {
      sz = (t->ch[0]) ? t->ch[0]->size : 0;
      if (k == sz + 1) break;
      if (k <= sz) t = t->ch[0];
      else k -= sz + 1, t = t->ch[1];
    }
    Splay(t, target);
  }
  // 前驱
  Node *predec(Node *x) {
    Node *u = x->ch[0];
    while (u && u->ch[1]) u = u->ch[1];
    return u;
  }
  // 后继
  Node *succes(Node *x) {
    Node *u = x->ch[1];
    while (u && u->ch[0]) u = u->ch[0];
    return u;
  }
}ST;