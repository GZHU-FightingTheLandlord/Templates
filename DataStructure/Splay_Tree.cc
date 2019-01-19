#include <bits/stdc++.h>
using namespace std;
 
const int maxn = 2e5 + 5;
 
struct Node {
  int key, rev, size;
  Node *fa, *ch[2];
  Node(int v = 0, Node *pre = NULL) : key(v), fa(pre) {
    ch[0] = ch[1] = NULL, rev = 0;
  }
  void pushup() {
    size = 1;
    if (ch[0] != NULL) size += ch[0]->size;
    if (ch[1] != NULL) size += ch[1]->size;
  }
  void pushdown() {
    if (rev) {
      if (ch[0]) ch[0]->reverse();
      if (ch[1]) ch[1]->reverse();
      rev = 0;
    }
  }
  void reverse() {
    swap(ch[0], ch[1]), rev ^= 1;
  }
};
 
struct SplayTree {
  int count;
  Node *root, box[maxn];
 
  void init() {
    root = NULL, count = 0;
  }
  Node *newnode(int key, Node *fa) {
    return &(box[++count] = Node(key, fa));
    // return new Node(key, fa);
  }
  // arr[0]与arr[n+1]需用作哨兵
  // build(ST.root, NULL, 0, n + 1, arr)
  void build(Node *&u, Node *fa, int l, int r, int *arr) {
    if (l > r) return (void)(u = NULL);
    int mid = (l + r) / 2;
    u = newnode(arr[mid], fa);
    build(u->ch[0], u, l, mid - 1, arr);
    build(u->ch[1], u, mid + 1, r, arr);
    u->pushup();
  }
  // 插入一个键值为v的结点，并旋转至树根，满足二叉搜索树性质
  void insert(int v) {
    if (root == NULL) {
      root = newnode(v, NULL);
    } else {
      Node *u = root; u->pushdown();
      while (u->ch[u->key < v]) {
        u->pushdown(), u = u->ch[u->key < v];
      }
      u->ch[u->key < v] = newnode(v, u);
      Splay(u->ch[u->key < v], NULL);
    }
  }
  // 删除结点x
  void erase(Node *x) {
    if (x == NULL) return;
    Splay(x, NULL);
    if (x->ch[0] && x->ch[1]) {
      Splay(succes(x), x);
      x->ch[1]->ch[0] = x->ch[0];
      x->ch[0]->fa = x->ch[1];
      x->ch[1]->fa = NULL, root = x->ch[1];
    } else if (x->ch[0]) {
      root = x->ch[0], x->ch[0]->fa = NULL;
    } else if (x->ch[1]) {
      root = x->ch[1], x->ch[1]->fa = NULL;
    } else {
      root = NULL;
    }
    // delete x;
  }
  // zig when c = 0, zag when c = 1
  void rotate(Node *x, int c) {
    Node *y = x->fa, *z = y->fa;
    y->pushdown(), x->pushdown();
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
    x->pushdown();
    while (x->fa != target) {
      Node *y = x->fa, *z = y->fa;
      if (z == target) {
        rotate(x, y->ch[0] == x);
      } else {
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
      t->pushdown();
      sz = (t->ch[0]) ? t->ch[0]->size : 0;
      if (k == sz + 1) break;
      if (k <= sz) t = t->ch[0];
      else k -= sz + 1, t = t->ch[1];
    }
    Splay(t, target);
  }
  // 前驱
  Node *predec(Node *x) {
    x->pushdown();
    Node *u = x->ch[0];
    while (u && u->ch[1]) {
      u->pushdown(), u = u->ch[1];
    }
    return u;
  }
  // 后继
  Node *succes(Node *x) {
    x->pushdown();
    Node *u = x->ch[1];
    while (u && u->ch[0]) {
      u->pushdown(), u = u->ch[0];
    }
    return u;
  }
  // 反转[l, r]
  void reverse(int l, int r) {
    if (l >= r) return;
    select(l, NULL);
    select(r + 2, root);
    Node *t = root->ch[1]->ch[0];
    t->reverse();
  }
  // 中序遍历
  void travel(Node *u, vector<int> &ans) {
    if (u) {
      u->pushdown();
      if (u->ch[0]) travel(u->ch[0], ans);
      if (u->key != -1) ans.push_back(u->key);
      if (u->ch[1]) travel(u->ch[1], ans);
    }
  }
  void print() {
    vector<int> ans;
    travel(root, ans);
    for (auto it = ans.begin(); it != ans.end(); it++) {
      if (it != ans.begin()) printf(" ");
      printf("%d", *it);
    }
    printf("\n");
  }
  // something else.
}ST;
