struct SAM {
  int last, tot, sz[maxn << 1], len[maxn << 1], fa[maxn << 1];
  int ch[maxn << 1][30];

  SAM() {
    tot = 0, last = newNode(0), len[0] = -1;
    memset(sz, 0, sizeof sz);
  }
  inline int newNode(int v) {
    len[++tot] = v, fa[tot] = 0;
    memset(ch[tot], 0, sizeof ch[tot]);
    return tot;
  }
  void append(int c) {
    int p = last, u = newNode(len[last] + 1);
    for (; p && !ch[p][c]; p = fa[p]) {
      ch[p][c] = u;
    }
    if (p == 0) {
      fa[u] = 1;
    } else {
      int q = ch[p][c];
      if (len[q] == len[p] + 1) {
        fa[u] = q;
      } else {
        int nq = newNode(len[p] + 1);
        memcpy(ch[nq], ch[q], sizeof ch[q]);
        fa[nq] = fa[q], fa[u] = fa[q] = nq;
        for (; p && (ch[p][c] == q); p = fa[p]) {
          ch[p][c] = nq;
        }
      }
    }
    last = u;
  }
  void match(char *s) {
    int pos = 1, length = 0;
    for (int i = 0, n = strlen(s); i < n; i++) {
      while (pos && !ch[pos][s[i] - 'a']) {
        pos = fa[pos], length = len[pos];
      }
      if (pos) {
        ++length, pos = ch[pos][s[i] - 'a'];
        // update ans
      } else {
        pos = 1, length = 0;
      }
    }
  }
} sam;