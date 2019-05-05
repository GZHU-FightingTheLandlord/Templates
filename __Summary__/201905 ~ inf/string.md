## String

### kmp

```cpp
template <template<class...> class T, class t>
VI getfail(const T<t>& s) {
  int n = sz(s);
  VI fail(n + 1);
  for (int i = 0, j = fail[0] = -1; i < n; i++, j++) {
    while (~j && s[j] != s[i]) j = fail[j];
    fail[i + 1] = j + 1;
  }
  return fail;
}

template <template<class...> class T, class t>
int match(const T<t> &s, const T<t> &par, const VI &fail) {
  int n = sz(s), m = sz(par);
  for (int i = 0, j = 0; i < n; ) {
    while (~j && par[j] != s[i]) j = fail[j];
    ++i, ++j;
    if (j >= m) return i - m + 1;
  }
  return -1;
}
```

### Z-function

```cpp
VI Zfunc(string s) {
  int n = sz(s);
  VI z(n);
  for (int i = 1, l = 0, r = 0; i < n; i++) {
    if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) ++z[i];
    if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
  }
  return z;
}
```

### Suffix Array

```cpp
template <typename T>
VI build_sa(int n, const T &s, int charset) {
  VI a(n);
  if (n == 0) {
    return a;
  }
  if (charset != -1) {
    VI aux(charset, 0);
    for (int i = 0; i < n; i++) {
      aux[ s[i] ]++;
    }
    int sum = 0;
    for (int i = 0; i < charset; i++) {
      int add = aux[i];
      aux[i] = sum;
      sum += add;
    }
    for (int i = 0; i < n; i++) {
      a[ aux[ s[i] ]++ ] = i;
    }
  } else {
    iota(a.begin(), a.end(), 0);
    sort(a.begin(), a.end(), [&s](int i, int j) { return s[i] < s[j]; });
  }
  VI sorted_by_second(n), ptr_group(n);
  VI new_group(n), group(n);
  group[ a[0] ] = 0;
  for (int i = 1; i < n; i++) {
    group[ a[i] ] = group[ a[i - 1] ] + ( !(s[ a[i] ] == s[ a[i - 1] ]) );
  }
  int cnt = group[a[n - 1]] + 1;
  int step = 1;
  while (cnt < n) {
    int at = 0;
    for (int i = n - step; i < n; i++) {
      sorted_by_second[at++] = i;
    }
    for (int i = 0; i < n; i++) {
      if (a[i] - step >= 0) {
        sorted_by_second[at++] = a[i] - step;
      }
    }
    for (int i = n - 1; i >= 0; i--) {
      ptr_group[ group[ a[i] ] ] = i;
    }
    for (int i = 0; i < n; i++) {
      int x = sorted_by_second[i];
      a[ ptr_group[ group[x] ]++ ] = x;
    }
    new_group[a[0]] = 0;
    for (int i = 1; i < n; i++) {
      if (group[ a[i] ] != group[ a[i - 1] ]) {
        new_group[ a[i] ] = new_group[ a[i - 1] ] + 1;
      } else {
        int pre = ( (a[i - 1] + step >= n) ? -1 : group[ a[i - 1] + step ] );
        int cur = ( (a[i] + step >= n) ? -1 : group[ a[i] + step ] );
        new_group[ a[i] ] = new_group[ a[i - 1] ] + (pre != cur);
      }
    }
    swap(group, new_group);
    cnt = group[ a[n - 1] ] + 1;
    step <<= 1;
  }
  return a;
}
```

one more

```cpp
namespace SuffixArray {
  const int maxn = "edit";

  int wa[maxn], wb[maxn], c[maxn], d[maxn];
  
  inline bool cmp(int *r, int a, int b, int k) {
    return (r[a] == r[b]) && (r[a + k] == r[b + k]);
  }
  
  void da(int *r, int *sa, int n, int m) {
    int i, j, p, *x = wa, *y = wb, *t;

    for (i = 0; i < m; i++) d[i] = 0;
    for (i = 0; i < n; i++) d[x[i] = r[i]]++;

    for (i = 1; i < m; i++) d[i] += d[i - 1];
    for (i = n - 1; i >= 0; i--) sa[--d[x[i]]] = i;

    for (j = 1, p = 1; j <= n; j <<= 1, m = p) {
      for (p = 0, i = n - j; i < n; i++) y[p++] = i;
      for (i = 0; i < n; i++) if (sa[i] >= j) y[p++] = sa[i] - j;

      for (i = 0; i < n; i++) c[i] = x[y[i]];
      for (i = 0; i < m; i++) d[i] = 0;

      for (i = 0; i < n; i++) d[c[i]]++;
      for (i = 1; i < m; i++) d[i] += d[i - 1];

      for (i = n - 1; i >= 0; i--) sa[--d[c[i]]] = y[i];
      for (t = x, x = y, y = t, p = 1, x[sa[0]] = 0, i = 1; i < n; i++) {
        x[sa[i]] = cmp(y, sa[i - 1], sa[i], j) ? (p - 1) : (p++);
      }
    }
  }
  int rank[maxn], height[maxn];
  void calheight(int *r, int *sa, int n) {
    int i, j, k = 0;
    for (i = 1; i <= n; i++) rank[sa[i]] = i;
    for (i = 0; i < n; i++) {
      if (k) --k;
      for (j = sa[rank[i] - 1]; r[i + k] == r[j + k]; k++);
      // blank
      height[rank[i]] = k;
    }
  }
}
```

### Suffix Automa

```cpp
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
```

### 最小表示法

对于一个字符串S，求S的循环的同构字符串S’中字典序最小的一个。

字符串"abcd"的循环同构字符串有：["abcd", "bcda", "cdab", "dabc"]。

```cpp
int minPresentation(string &s) {
  int n = s.length();
  int i = 0, j = 1, k = 0;
  while (k < n && i < n && j < n) {
    if (s[(i + k) % n] == s[(j + k) % n]) {
      ++k;
    } else {
      s[(i + k) % n] > s[(j + k) % n] ? (i += k + 1) : (j += k + 1);
      i += (i == j);
      k = 0;
    }
  }
  return min(i, j);
}
```

### Manacher

```cpp
int p[maxn << 1];
char str[maxn << 1];

int manacher(char *s, int n) {
  str[0] = '$'; str[1] = '#';

  for (int i = 0; i < n; i++) {
    str[(i << 1) + 2] = s[i];
    str[(i << 1) + 3] = '#';
  }
  n = (n + 1) << 1;
  str[n] = 0;

  int ret = 0, mx = 0, pos;
  for (int i = 1; i < n; i++) {
    p[i] = mx > i ? min(p[(pos << 1) - i], mx - i) : 1;

    while (str[i - p[i]] == str[i + p[i]]) p[i]++;

    if (p[i] + i > mx) mx = p[i] + i, pos = i;

    ret = max(ret, p[i]);
  }
  return ret - 1;
}
```

### AC Automa

```cpp
namespace acam {
  struct Node {
    int son[26], fail;
    void init() {
      fail = 0;
      memset(son, 0, sizeof son);
    }
  } T[N];
  int tot;

#define Son(i, x) T[i].son[x]
#define trans(c) (c - 'A')

  void init() {
    tot = 0, T[0].init();
  }
  void insert(char *s, int index) {
    int cur = 0;
    for (int i = 0; s[i]; i++) {
      int c = trans(s[i]);
      if (!Son(cur, c)) {
        Son(cur, c) = ++tot;
        T[tot].init();
      }
      cur = Son(cur, c);
    }
  }
  void build() {
    queue<int> Q;
    for (int i = 0; i < 26; i++) {
      if (Son(0, i)) Q.push(Son(0, i));
    }
    while (!Q.empty()) {
      int u = Q.front(); Q.pop();
      for (int i = 0; i < 26; i++) {
        if (Son(u, i)) {
          T[Son(u, i)].fail = T[T[u].fail].son[i];
          Q.push(Son(u, i));
        } else {
          T[u].son[i] = T[T[u].fail].son[i];
        }
      }
    }
  }
  int query(char *t) {
    int ans = 0, cur = 0;
    for (int i = 0; t[i]; i++) {
      int c = trans(t[i]);
      cur = Son(cur, c);
      for (int j = cur; j; j = T[j].fail) {
        // upd ans;
      }
    }
    return ans;
  }
}
```

### PAM

```cpp
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
```