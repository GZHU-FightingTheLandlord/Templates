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
// sais O(n)
// unfinished but usable

// 可以对solve传string或vector<int>等
// 用完solve之后用sa rk ht即可
namespace SA {
  const size_t sz = 3e5 + 5;
  int bucket[sz], bucket1[sz], sa[sz], rk[sz], ht[sz];
  bool type[sz << 1];
  bool isLMS(const int i, const bool *type) {
    return i > 0 && type[i] && !type[i - 1];
  }

  template<class T>
  void inducedSort(const T &s, int *sa, const int len, const int sigma, const int bucketSize, bool *type, int *bucket, int *cntbuf, int *p) {
    memset(bucket, 0, sizeof(int) * sigma);
    memset(sa, -1, sizeof(int) * len);
    for (int i = 0; i < len; i++) {
      bucket[s[i]]++;
    }
    cntbuf[0] = bucket[0];
    for (int i = 1; i < sigma; i++) {
      cntbuf[i] = cntbuf[i - 1] + bucket[i];
    }
    for (int i = bucketSize - 1; i >= 0; i--) {
      sa[--cntbuf[s[p[i]]]] = p[i];
    }
    for (int i = 1; i < sigma; i++) {
      cntbuf[i] = cntbuf[i - 1] + bucket[i - 1];
    }
    for (int i = 0; i < len; i++) {
      if (sa[i] > 0 && !type[sa[i] - 1]) {
        sa[cntbuf[s[sa[i] - 1]]++] = sa[i] - 1;
      }
    }
    cntbuf[0] = bucket[0];
    for (int i = 1; i < sigma; i++) {
      cntbuf[i] = cntbuf[i - 1] + bucket[i];
    }
    for (int i = len - 1; i >= 0; i--) {
      if (sa[i] > 0 && type[sa[i] - 1]) {
        sa[--cntbuf[s[sa[i] - 1]]] = sa[i] - 1;
      }
    }
  }
  template<class T>
  void sais(const T &s, int *sa, int len, bool *type, int *bucket, int *bucket1, int sigma) {
    int i, j, bucketSize = 0, cnt = 0, p = -1, x, *cntbuf = bucket + sigma;
    type[len - 1] = 1;
    for (i = len - 2; i >= 0; i--) {
      type[i] = s[i] < s[i + 1] || (s[i] == s[i + 1] && type[i + 1]);
    }
    for (i = 1; i < len; i++) {
      if (type[i] && !type[i - 1]) {
        bucket1[bucketSize++] = i;
      }
    }
    inducedSort(s, sa, len, sigma, bucketSize, type, bucket, cntbuf, bucket1);
    for (i = bucketSize = 0; i < len; i++) {
      if (isLMS(sa[i], type)) {
        sa[bucketSize++] = sa[i];
      }
    }
    for (i = bucketSize; i < len; i++) {
      sa[i] = -1;
    }
    for (i = 0; i < bucketSize; i++) {
      x = sa[i];
      for (j = 0; j < len; j++) {
        if (p == -1 || s[x + j] != s[p + j] || type[x + j] != type[p + j]) {
          cnt++, p = x;
          break;
        } else if (j > 0 && (isLMS(x + j, type) || isLMS(p + j, type))) {
          break;
        }
      }
      x = (~x & 1 ? x >> 1 : (x - 1) >> 1), sa[bucketSize + x] = cnt - 1;
    }
    for (i = j = len - 1; i >= bucketSize; i--) {
      if (sa[i] >= 0) {
        sa[j--] = sa[i];
      }
    }
    int *s1 = sa + len - bucketSize, *bucket2 = bucket1 + bucketSize;
    if (cnt < bucketSize) {
      sais(s1, sa, bucketSize, type + len, bucket, bucket1 + bucketSize, cnt);
    } else {
      for (i = 0; i < bucketSize; i++) {
        sa[s1[i]] = i;
      }
    }
    for (i = 0; i < bucketSize; i++) {
      bucket2[i] = bucket1[sa[i]];
    }
    inducedSort(s, sa, len, sigma, bucketSize, type, bucket, cntbuf, bucket2);
  }

  void getHeight(const vector<int> &s, int n) {
    int i, j, k = 0;
    for(i = 1; i <= n; i++) {
      rk[sa[i]] = i;
    }
    for(i = 0; i < n; i++) {
      if(k) k--;
      for(j = sa[rk[i] - 1]; s[i + k] == s[j + k]; k++);
      ht[rk[i]] = k;
    }
  }

  template<class T>
  void solve(const T &s) {
    const int n = s.size();
    vector<int> v(n + 1);
    int sigma = 0;
    for(int i = 0; i < n; i++) {
      v[i] = int(s[i]) + 1;
      sigma = max(sigma, v[i] + 1);
    }
    v[n] = 0;
    sais(v, sa, n + 1, type, bucket, bucket1, sigma);
    getHeight(v, n);
  }
}