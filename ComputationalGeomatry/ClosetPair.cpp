struct ClosetPair {
  static const int N = 1e5 + 10;
  struct pt {
    double x, y;
    int id;
  };
  vector<pt> a, t;
  void init() {
    a.clear();
  }
  static int dcmp(double x, double y) {
    if(fabs(x - y) < 1e-7) return 0;
    return x < y ? -1 : 1;
  }
  void addPoint(double x, double y) {
    pt psh;
    tie(psh.x, psh.y, psh.id) = make_tuple(x, y, int(a.size()));
    a.push_back(psh);
  }
  double minDist;
  int a1 = -1, a2 = -1;
  double sq(double x) {
    return x * x;
  }
  void upd(const pt& x, const pt& y) {
    double dis = sqrt(sq(x.x - y.x) + sq(x.y - y.y));
    if(dis < minDist) {
      tie(minDist, a1, a2) = make_tuple(dis, x.id, y.id);
    }
  }
  void _solve(int l, int r) {
    if(r - l <= 3) {
      for(int i = l; i <= r; i++) {
        for(int j = i + 1; j <= r; j++) {
          upd(a[i], a[j]);
        }
      }
      sort(&a[l], &a[r + 1], [&](const pt& x, const pt& y) {
        return dcmp(x.y, y.y) < 0;
      });
      return;
    }
    int m = (l + r) >> 1;
    double midx = a[m].x;
    _solve(l, m);
    _solve(m + 1, r);
    merge(&a[l], &a[m + 1], &a[m + 1], &a[r + 1], &t[0], [&](const pt& x, const pt& y) {
      return dcmp(x.y, y.y) < 0;
    });
    copy(&t[0], &t[r - l + 1], &a[l]);
    int tsz = 0;
    for(int i = l; i <= r; i++) {
      if(abs(a[i].x - midx) < minDist) {
        for(int j = tsz - 1; j >= 0 && a[i].y - t[j].y < minDist; j--) {
          upd(a[i], t[j]);
        }
        t[tsz++] = a[i];
      }
    }
  }
  void solve() {
    sort(a.begin(), a.end(), [&](const pt& x, const pt& y) {
      int cmp = dcmp(x.x, y.x);
      return cmp < 0 || (cmp == 0 && dcmp(x.y, y.y) < 0);
    });
    tie(minDist, a1, a2) = make_tuple(1e200, -1, -1);
    const int n = a.size();
    t.resize(n);
    _solve(0, n - 1);
  }
};
