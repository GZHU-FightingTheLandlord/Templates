#include <bits/stdc++.h>
using namespace std;

struct leftist_heap {
    static const int maxn = 2e5 + 5;
    static const int maxm = 1e3 + 5;

    int root[maxm], tot;
    int v[maxn], ch[maxn][2], size[maxn], stk[maxn], tp;

    inline void init(int &n) {
        tot = tp = 0;
        memset(root, 0, (n + 1) * sizeof(int));
    }

    inline int newnode(int x) {
        int rt = (tp > 0) ? stk[--tp] : (tot++);
        v[rt] = x, size[rt] = 1, ch[rt][0] = ch[rt][1] = 0;
        return rt;
    }

    int merge(int a, int b) {
        if (!a || !b) {
            return a ? a : b;
        }
        if (v[a] > v[b]) {
            swap(a, b);
        }
        ch[a][1] = merge(ch[a][1], b);
        if (size[ch[a][0]] < size[ch[a][1]]) {
            swap(ch[a][0], ch[a][1]);
        }
        size[a] = size[ch[a][1]] + 1;
        return a;
    }

    void Merge(int a, int b) {
        root[a] = merge(root[a], root[b]);
        root[b] = 0;
    }

    void insert(int i, int val) {
        int tmp = newnode(val);
        root[i] = merge(root[i], tmp);
    }

    inline int top(int i) {
        return (root[i] > 0) ? v[root[i]] : -1;
    }

    inline void pop(int i) {
        (root[i] > 0) ? (stk[tp++] = root[i], root[i] = merge(ch[root[i]][0], ch[root[i]][1])) : 0;
    }
};