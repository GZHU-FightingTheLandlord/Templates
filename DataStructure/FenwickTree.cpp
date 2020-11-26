#include <bits/stdc++.h>

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


#if 0

namespace data_structure {

template<typename Type = int>
class FenwickTree {
 public:
  explict FenwickTree(const int size) {
    n_ = size;
    fenwick_tree_.resize(size + 1);
  }
  explict FenwickTree(const int size, const Type& init_value) {
    std::vector<Type> vec(size + 1, init_value);
    Init(vec);
  }
  explicit FenwickTree(const vector<Type>& vec) {
    Init(vec);
  }
  void Add(int position, const Type& value) {
    assert(position > 0);
    for (int i = position; i <= n_; i += i & (-i)) {
      fenwick_tree_[i] += value;
    }
  }
  Type Query(int index) const {
    Type ret;
    for (int i = index; i > 0; i -= i & (-i)) {
      ret += fenwick_tree_[i];
    }
    return ret;
  }
  Type Query(int left, int right) const {
    assert(left > 0 && right >= left);
    return Query(right) - Query(left - 1);
  }
 private:
  void Init(const std::vector<Type>& vec) {
    n_ = vec.size() - 1;
    fenwick_tree_.resize(vec.size());
    for (int i = 1; i <= n_; i++) {
      Add(i, vec[i]);
    }
  }
  int n_;
  std::vector<Type> fenwick_tree_;
};

}

#endif