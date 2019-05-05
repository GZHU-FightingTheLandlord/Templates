struct PAM {
  struct Node {
    int son[27], fail, len, dep;
    void init(int l, int f, int d = 0) {
      fail = f, len = l, dep = d;
      memset(son, 0, sizeof son);
    }
  } T[N + 5];
  int tot, prefix, suffix, l, r;
  int s[N * 2 + 5];

  void init() {
    ans = 0;
    tot = 1, l = N + 1, r = N, prefix = suffix = 0;
    T[0].init(0, 1), T[1].init(-1, 0);
    memset(s, 0, sizeof s);
  }
  void encode(int &c) {
    // keep c > 0
    c = c - 'a' + 1;
  }
  int pre_fail(int cur) {
    while (s[l + T[cur].len + 1] != s[l]) {
      cur = T[cur].fail;
    }
    return cur;
  }
  int suf_fail(int cur) {
    while (s[r - T[cur].len - 1] != s[r]) {
      cur = T[cur].fail;
    }
    return cur;
  }
  void push_front(int c) {
    encode(c), s[--l] = c;
    prefix = pre_fail(prefix);
    if (!T[prefix].son[c]) {
      int f = pre_fail(T[prefix].fail);
      T[++tot].init(T[prefix].len + 2, T[f].son[c], T[T[f].son[c]].dep + 1);
      T[prefix].son[c] = tot;
    }
    prefix = T[prefix].son[c];
    if (T[prefix].len == r - l + 1) {
      suffix = prefix;
    }
  }
  void push_back(int c) {
    encode(c), s[++r] = c;
    suffix = suf_fail(suffix);
    if (!T[suffix].son[c]) {
      int f = suf_fail(T[suffix].fail);
      T[++tot].init(T[suffix].len + 2, T[f].son[c], T[T[f].son[c]].dep + 1);
      T[suffix].son[c] = tot;
    }
    suffix = T[suffix].son[c];
    if (T[suffix].len == r - l + 1) {
      prefix = suffix;
    }
  }
} pam;