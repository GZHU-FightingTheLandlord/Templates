struct Seive {
  int maxn;
  vector<bool> isp;
  vector<int> p, phi, mu;

  Seive(int n = 0) : maxn(n), isp(n + 5, true), phi(n + 5, 0), mu(n + 5, 0) { solve(); }

  void solve() {
    isp[0] = isp[1] = false; phi[1] = 1; mu[1] = 1;
    for (int i = 2; i <= maxn; i++) {
      if (isp[i]) {
        p.push_back(i);
        phi[i] = i - 1;
        mu[i] = -1;
      }
      for (int j = 0; j < (int)p.size() && i * p[j] <= maxn; j++) {
        const int cur = i * p[j];
        isp[cur] = false;
        if (i % p[j]) {
          phi[cur] = phi[i] * (p[j] - 1);
          mu[cur] = -mu[i];
        }
        else {
          phi[cur] = phi[i] * p[j];
          mu[cur] = 0;
          break;
        }
      }
    }
  }
};
