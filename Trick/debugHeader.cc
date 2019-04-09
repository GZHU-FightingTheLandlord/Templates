#ifdef yuki
#define debug(x...) do { cout << "["#x"] -> "; err(x); } while (0)
template <class T>
inline void _E(const T& x) { cout << x; }
inline void _E(const double& x) {
  string str = to_string(x);
  cout << str.substr(0, str.find('.')) + str.substr(str.find('.'), 4);
}
inline void _E(const long double& x) { _E((double)x); }
template <class T>
inline void _E(const vector<T>& vec) {
  for (auto it = begin(vec), en = end(vec); it != en; it++) {
    cout << (it != begin(vec) ? " " : ""); _E(*it);
  }
}
inline void err() { cout.flush(); }
template <class T, class... U>
inline void err(const T& arg, const U&... args) {
  cout << "["; _E(arg); cout << (sizeof...(args) ? "], " : "]\n"); err(args...);
}
#else
#define debug(...) do {} while (0)
#endif