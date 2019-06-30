// 能用但未整理 考完试再来

namespace SA {
  #define F(x) ((x)/3+((x)%3==1?0:tb))
  #define G(x) ((x)<tb?(x)*3+1:((x)-tb)*3+2)
  #define rep(i, a, b) for(int i = a; i < b; i++)
  #define per(i, a, b) for(int i = b - 1; i >= a; i--)
  const int N = 2e4 + 5;
  const int MAXN = N * 10;
  int wv[MAXN], wws[MAXN], wa[MAXN], wb[MAXN], sa[MAXN], r[MAXN], height[MAXN], rank[MAXN];
  void Sort(int *r, int *a, int *b, int n, int m) {
    memset(wws, 0, sizeof(int) * m);
    rep(i, 0, n) wv[i] = r[a[i]];
    rep(i, 0, n) wws[wv[i]]++;
    rep(i, 1, m) wws[i] += wws[i - 1];
    per(i, 0, n) b[--wws[wv[i]]] = a[i];
  }
  int c0(int *r, int a, int b) {
    return r[a] == r[b] && r[a + 1] == r[b + 1] && r[a + 2] == r[b + 2];
  }
  int c12(int k, int *r, int a, int b) {
    if(k == 2) return r[a] < r[b] || r[a] == r[b] && c12(1, r, a + 1, b + 1);
    return r[a] < r[b] || r[a] == r[b] && wv[a + 1] < wv[b + 1];
  }
  void dc3(int *r, int *sa, int n, int m) {
    n++, m++;
    int i, j, *rn = r + n, *san = sa + n, ta = 0, tb = (n + 1) / 3, tbc = 0, p;
    r[n] = r[n + 1] = 0;
    rep(i, 0, n) if(i % 3 != 0) wa[tbc++] = i;
    Sort(r + 2, wa, wb, tbc, m);
    Sort(r + 1, wb, wa, tbc, m);
    Sort(r, wa, wb, tbc, m);
    for(p = 1, rn[F(wb[0])] = 0, i = 1; i < tbc; i++) rn[F(wb[i])] = c0(r, wb[i - 1], wb[i]) ? p - 1 : p++;
    if(p < tbc) dc3(rn, san, tbc - 1, p - 1);
    else rep(i, 0, tbc) san[rn[i]] = i;
    rep(i, 0, tbc) if(san[i] < tb) wb[ta++] = san[i] * 3;
    if(n % 3 == 1) wb[ta++] = n - 1;
    Sort(r, wb, wa, ta, m);
    rep(i, 0, tbc) wv[wb[i] = G(san[i])] = i;
    for(i = 0, j = 0, p = 0; i < ta && j < tbc; p++) sa[p] = c12(wb[j] % 3, r, wa[i], wb[j]) ? wa[i++] : wb[j++];
    for(; i < ta; p++) sa[p] = wa[i++];
    for(; j < tbc; p++) sa[p] = wb[j++];
  }
  void calheight(int *r, int *sa, int n) {
    int k = 0;
    rep(i, 1, n + 1) rank[sa[i]] = i;
    for(int i = 0; i < n; height[rank[i++]] = k) {
      if(k) k--;
      for(int j = sa[rank[i] - 1]; r[i + k] == r[j + k]; k++);
    }
  }
  // dc3(SA::r, SA::sa, SA::r.size(), max{SA::r});
  // calheight(SA::r, SA::sa, SA::r.size());
};
