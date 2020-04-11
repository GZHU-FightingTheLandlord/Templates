class strHash {
  vector<pair<int, int>> modAndSeed; // {mod, seed}
  vector<vector<int>> invPow, pre;
  int sz;
  string str;
  int fpow(int a, int b, int c) {
    int ans = 1;
    for(a %= c; b > 0; b >>= 1, a = 1ll * a * a % c) {
      if(b & 1) ans = 1ll * ans * a % c;
    }
    return ans;
  }
  void add(int &x, int y, int mod) {
    x += y;
    if(x >= mod) x -= mod;
  }
  void sub(int &x, int y, int mod) {
    x -= y;
    if(x < 0) x += mod;
  }
  void insertModAndSeed() {}
  template<typename T, typename... U>
  void insertModAndSeed(const T& arg, const U&... args) {
    modAndSeed.push_back(arg);
    const int inv = fpow(arg.second, arg.first - 2, arg.first);
    invPow.push_back(vector<int>(sz + 1, 1));
    for(int i = 1; i <= sz; i++) {
      invPow.back()[i] = 1ll * invPow.back()[i - 1] * inv % arg.first;
    }
    pre.push_back(vector<int>(sz + 1, 0));
    int curSeed = arg.second;
    for(int i = 0; i < sz; i++) {
      pre.back()[i + 1] = pre.back()[i];
      add(pre.back()[i + 1], 1ll * curSeed * int(str[i]) % arg.first, arg.first);
      curSeed = 1ll * arg.second * curSeed % arg.first;
    }
    insertModAndSeed(args...);
  }
  void defaultModAndSeed() {
    mt19937 obj(chrono::system_clock::now().time_since_epoch().count());
    const int m1 = 998434903, m2 = 1019991001;
    const int s1 = obj() >> 13, s2 = obj() >> 13;
    insertModAndSeed(make_pair(m1, s1), make_pair(m2, s2));
  }
  int get(int id, int l, int r) {
    int ret = pre[id][r];
    sub(ret, pre[id][l - 1], modAndSeed[id].first);
    ret = 1ll * ret * invPow[id][l - 1] % modAndSeed[id].first;
    return ret;
  }
public:
  struct node {
    vector<int> a;
    int& operator[](const int& x) { return a[x]; }
    const int& operator[](const int& x) const { return a[x]; }
    bool operator<(const node& rhs) const {
      for(int i = 0; i < int(a.size()); i++) {
        if(a[i] < rhs[i]) return true;
        if(a[i] > rhs[i]) return false;
      }
      return false;
    }
    bool operator==(const node& rhs) const {
      for(int i = 0; i < int(a.size()); i++) {
        if(a[i] != rhs[i]) return false;
      }
      return true;
    }
  };
  template<typename... U>
  strHash(const string& s, const U&... args) {
    str = s;
    sz = s.size();
    insertModAndSeed(args...);
    if(modAndSeed.empty()) defaultModAndSeed();
  }
  node get(int l, int r) {
    l++, r++;
    node ret;
    const int n = modAndSeed.size();
    ret.a.resize(n);
    for(int i = 0; i < n; i++) {
      ret[i] = get(i, l, r);
    }
    return ret;
  }
};
