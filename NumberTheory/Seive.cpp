struct Seive {
  int n;
  vector<bool> isp;
  vector<int> p, phi, mu, sig, num;
  Seive(int nn) : n(nn), isp(nn, true), phi(nn, 0), mu(nn, 0), sig(nn, 0), num(nn, 0) { solve(); }
  void solve() {
    isp[0] = isp[1] = false; phi[1] = 1; mu[1] = 1; sig[1] = 1;
    for (int i = 2; i < n; i++) {
      if (isp[i]) {
        p.push_back(i);
        phi[i] = i - 1;
        mu[i] = -1;
        num[i] = 1;
        sig[i] = 2;
      }
      for (int j = 0; j < (int)p.size() && i * p[j] < n; j++) {
        const int cur = i * p[j];
        isp[cur] = false;
        if (i % p[j]) {
          num[cur] = 1;
          sig[cur] = sig[i] * sig[p[j]];
          phi[cur] = phi[i] * (p[j] - 1);
          mu[cur] = -mu[i];
        }
        else {
          num[cur] = num[i] + 1;
          sig[cur] = sig[i] / (num[i] + 1) * (num[i] + 2);
          phi[cur] = phi[i] * p[j];
          mu[cur] = 0;
          break;
        }
      }
    }
  }
};
