// ~i & i + 1
// i & -i

struct bit {
  int N, base[maxn];
  void setN(int n) { N = n; }
  void init() { memset(base, 0, sizeof base); }
  void add(int at, int v) {
    if (!at) return;
    for (int i = at; i <= N; i += i & -i) {
      base[i] += v;
    }
  }
  int getSum(int at) {
    int sum = 0;
    for (int i = at; i; i -= i & -i) {
      sum += base[i];
    }
    return sum;
  }
};

struct Interval {
  int N, base[2][maxn];
  void setN(int n) { N = n; }
  void init() { memset(base, 0, sizeof base); }
  void add(int at, int v) {
    if (!at) return;
    for (int i = at; i <= N; i += i & -i) {
      base[0][i] += v, base[1][i] -= v * at;
    }
  }
  void add(int l, int r, int v) {
    add(l, v), add(r + 1, -v);
  }
  int getSum(int at) {
    int sum = 0, mul = at + 1;
    for (int i = at; i; i -= i & -i) {
      sum += mul * base[0][i] + base[1][i];
    }
    return sum;
  }
  int query(int l, int r) {
    return getSum(r) - getSum(l - 1);
  }
};