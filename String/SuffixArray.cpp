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



// ****************************************************************************************
// dc3 O(n)
// unfinished but usable

namespace SA {
  #define F(x) ((x)/3+((x)%3==1?0:tb))
  #define G(x) ((x)<tb?(x)*3+1:((x)-tb)*3+2)
  const int N = 2e4 + 5;
  const int MAXN = N * 4; // maybe 3 times is enough?
  int wv[N], wws[N], wa[N], wb[N], ssa[MAXN], r[MAXN];
  void Sort(int *r, int *a, int *b, int n, int m) {
    memset(wws, 0, sizeof(int) * m);
    for(int i = 0; i < n; i++) wv[i] = r[a[i]];
    for(int i = 0; i < n; i++) wws[wv[i]]++;
    for(int i = 1; i < m; i++) wws[i] += wws[i - 1];
    for(int i = n - 1; i >= 0; i--) b[--wws[wv[i]]] = a[i];
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
    for(i = 0; i < n; i++) if(i % 3 != 0) wa[tbc++] = i;
    Sort(r + 2, wa, wb, tbc, m);
    Sort(r + 1, wb, wa, tbc, m);
    Sort(r, wa, wb, tbc, m);
    for(p = 1, rn[F(wb[0])] = 0, i = 1; i < tbc; i++) rn[F(wb[i])] = c0(r, wb[i - 1], wb[i]) ? p - 1 : p++;
    if(p < tbc) dc3(rn, san, tbc - 1, p - 1);
    else for(i = 0; i < tbc; i++) san[rn[i]] = i;
    for(i = 0; i < tbc; i++) if(san[i] < tb) wb[ta++] = san[i] * 3;
    if(n % 3 == 1) wb[ta++] = n - 1;
    Sort(r, wb, wa, ta, m);
    for(i = 0; i < tbc; i++) wv[wb[i] = G(san[i])] = i;
    for(i = 0, j = 0, p = 0; i < ta && j < tbc; p++) sa[p] = c12(wb[j] % 3, r, wa[i], wb[j]) ? wa[i++] : wb[j++];
    for(; i < ta; p++) sa[p] = wa[i++];
    for(; j < tbc; p++) sa[p] = wb[j++];
  }
  // below for whom don't know how to use the above
  template<typename T>
  vector<int> SA(const T &s) {
    const int n = s.size();
    int mx = 0;
    for(int i = 0; i < n; i++) {
      r[i] = int(s[i]) + 1;
      mx = max(r[i], mx);
    }
    dc3(r, ssa, n, mx);
    vector<int> ret(n + 1);
    for(int i = 0; i <= n; i++) {
      ret[i] = ssa[i];
    }
    return ret;
  }
  vector<int> RANK(const vector<int> &sa) {
    vector<int> ret(sa.size() - 1);
    for(int i = 1; i < int(sa.size()); i++) {
      ret[sa[i]] = i;
    }
    return ret;
  }
  template<typename T>
  vector<int> HEIGHT(const T &s, const vector<int> &sa, const vector<int> &rank) {
    const int n = int(rank.size());
    int k = 0; vector<int> ret(n + 1);
    for(int i = 0; i < n; ret[rank[i++]] = k) {
      if(k) k--;
      for(int j = sa[rank[i] - 1]; s[i + k] == s[j + k]; k++);
    }
    return ret;
  }
};
