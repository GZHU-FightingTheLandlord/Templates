// The part of calculate sa[] is fine, but height[] may not.

// sa[i]: if suffix(k) is the i-th smallest suffix, then sa[i] = k
// height[i]: lcp of suffix(sa[i]) and suffix(sa[i + 1])
template <int charset = 256> struct SA {
  vector<int> sa, height;
#define getrank(x) (((x) < (int)rank.size()) ? rank[(x)] : -1)
  SA(string s = "") {
    int n = s.size();
    sa.resize(n);
    vector<int> rank(n), cand(charset, 0);
    for (int i = 0; i < n; i++) {
      cand[s[i]]++;
    }
    for (int i = 1; i < charset; i++) {
      cand[i] += cand[i - 1];
    }
    for (int i = n - 1; i >= 0; i--) {
      sa[--cand[s[i]]] = i;
    }

    rank[sa[0]] = 0;
    for (int i = 1; i < n; i++) {
      rank[sa[i]] = rank[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
    }
    cand.resize(n + 5, 0);
    for (int l = 1; l < n; l <<= 1) {
      vector<int> newsa(n);
      fill(cand.begin(), cand.end(), 0);
      for (int i = 0; i < n; i++) {
        cand[getrank(sa[i] + l) + 1]++;
      }
      for (int i = 1; i < n + 5; i++) {
        cand[i] += cand[i - 1];
      }
      for (int i = n - 1; i >= 0; i--) {
        newsa[--cand[getrank(sa[i] + l) + 1]] = sa[i];
      }
      fill(cand.begin(), cand.end(), 0);
      for (int i = 0; i < n; i++) {
        cand[rank[newsa[i]]]++;
      }
      for (int i = 1; i < n + 5; i++) {
        cand[i] += cand[i - 1];
      }
      for (int i = n - 1; i >= 0; i--) {
        sa[--cand[rank[newsa[i]]]] = newsa[i];
      }

      vector<int> newrank(n, 0);
      for (int i = 1; i < n; i++) {
        newrank[sa[i]] = newrank[sa[i - 1]];
        if (rank[sa[i - 1]] != rank[sa[i]] || getrank(sa[i - 1] + l) != getrank(sa[i] + l)) {
          ++newrank[sa[i]];
        }
      }
      swap(rank, newrank);
    }

    for (int i = 0; i < n; i++) rank[sa[i]] = i;
    
    int k = 0;
    height.resize(max(n - 1, 0));
    for (int i = 0; i < n; i++) {
      k = max(k - 1, 0);
      if (rank[i] == n - 1) {
        k = 0;
      } else {
        int j = sa[rank[i] + 1];
        while (i + k < n && j + k < n && s[i + k] == s[j + k]) {
          ++k;
        }
        height[rank[i]] = k;
      }
    }
  }
#undef getrank
};