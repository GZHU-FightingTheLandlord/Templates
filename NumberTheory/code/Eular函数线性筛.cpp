#include <algorithm>
#include <vector>
using namespace std;

struct Seive {
	int maxn;
	vector<bool> isp;
	vector<int> p, phi;

	Seive(int n = 0) : maxn(n), isp(n + 5), phi(n + 5) { 
		for (int i = 0; i < n + 5; i++) {
			isp[i] = true;
			phi[i] = 0;
		}
		solve();
	}

	void solve() {
		isp[0] = isp[1] = false; phi[1] = 1;
		for (int i = 2; i <= maxn; i++) {
			if (isp[i]) p.push_back(i), phi[i] = i - 1;
			for (int j = 0; j < (int)p.size() && i * p[j] <= maxn; j++) {
				isp[i * p[j]] = false;
				if (i % p[j]) {
					phi[i * p[j]] = phi[i] * (p[j] - 1);
				}
				else {
					phi[i * p[j]] = phi[i] * p[j];
					break;
				}
			}
		}
	}

	inline bool chkpri(int x) { return isp[x]; }
	inline int getphi(int x) { return phi[x]; }
};