namespace GCD {
  const int N = 1000000; // O(N)
  const int SN = 1000; // SN = sqrt(N)
  bool np[N + 5];
  int ps[N + 5], cs[N + 5][3], pn = 0, _gcd[SN + 3][SN + 3];
  void init() {
    np[1] = true;
    cs[1][0] = cs[1][1] = cs[1][2] = 1;
    for(int i = 2; i <= N; i++) {
      if(!np[i]) {
        cs[i][0] = cs[i][1] = 1;
        cs[i][2] = i;
        ps[++pn] = i;
      }
      for(int j = 1; j <= pn && i * ps[j] <= N; j++) {
        np[i * ps[j]] = 1;
        int cm = cs[i][0] * ps[j];
        if(cm < cs[i][1]) {
          cs[i * ps[j]][0] = cm;
          cs[i * ps[j]][1] = cs[i][1];
          cs[i * ps[j]][2] = cs[i][2];
        } else if(cm < cs[i][2]) {
          cs[i * ps[j]][0] = cs[i][1];
          cs[i * ps[j]][1] = cm;
          cs[i * ps[j]][2] = cs[i][2];
        } else {
          cs[i * ps[j]][0] = cs[i][1];
          cs[i * ps[j]][1] = cs[i][2];
          cs[i * ps[j]][2] = cm;
        }
        if(i % ps[j] == 0) {
          break;
        }
      }
    }
    for(int i = 0; i <= SN; i++) {
      _gcd[i][0] = _gcd[0][i] = i;
    }
    for(int i = 1; i <= SN; i++) {
      for(int j = 1; j <= i; j++) {
        _gcd[i][j] = _gcd[j][i] = _gcd[i - j][j];
      }
    }
  }
  int gcd(int a, int b) {
    // assert(a <= N && b <= N);
    int *x = cs[a], g = 1;
    for(int i = 0; i < 3; i++) {
      int d;
      if(x[i] <= SN) d = _gcd[x[i]][b % x[i]];
      else if(b % x[i]) d = 1;
      else d = x[i];
      g *= d; b /= d;
    }
    return g;
  }
  static int __GCD_INIT = []() {
    init();
    return 0;
  }();
}
