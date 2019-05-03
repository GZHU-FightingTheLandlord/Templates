#define sz(x) ((int)x.size())
typedef vector<int> VI;

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