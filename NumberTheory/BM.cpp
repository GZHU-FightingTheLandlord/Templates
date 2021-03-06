/*使用说明：
输入项数n，和数列前几项，输出其最可能的线性递推方程

例子1：
输入：
7
1 1 2 3 5 8 13

输出：
1.000000 1.000000
该输出表示，递推式为f(x)=1*f(x-1)+1*f(x-2)

例子2：
输入：
9
1 2 3 4 5 6 7 8 9

输出：
2.000000 -1.000000
该输出表示递推式为f(x)=2*f(x-1)-1*f(x-2)

例子3：
输入：
9
1 1 1 2 2 3 3 4 4

输出：
1.000000 1.000000 -1.000000 0.000000
该输出表示递推式为f(x)=1*f(x-1)+1*f(x-2)-1*f(x-3)+0*f(x-4)
*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;
struct BM {
	static const int MAXN = 10005;
	int n, pn, fail[MAXN];
	double delta[MAXN];
	vector<double> ps[MAXN];
	void Solve(double x[], const int &n) {
		pn = 0;
		memset(fail, 0, sizeof fail);
		memset(delta, 0, sizeof delta);
		ps[0].clear();
		for (int i = 1; i <= n; i++) {
			double dt = -x[i];
			for (int j = 0; j < ps[pn].size(); j++) {
				dt += x[i - j - 1] * ps[pn][j];
			}
			delta[i] = dt;
			if (fabs(dt) <= 1e-8) continue;
			fail[pn] = i;
			if (!pn) {
				ps[++pn].resize(1);
				continue;
			}
			vector<double> &ls = ps[pn - 1];
			double k = -dt / delta[fail[pn - 1]];
			vector<double> cur(i - fail[pn - 1] - 1);
			cur.push_back(-k);
			for (int j = 0; j < ls.size(); j++) {
				cur.push_back(ls[j] * k);
			}
			if (cur.size() < ps[pn].size()) {
				cur.resize(ps[pn].size());
			}
			for (int j = 0; j < ps[pn].size(); j++) {
				cur[j] += ps[pn][j];
			}
			ps[++pn] = cur;
		}
	}
	void print() {
		for (int i = 0; i < ps[pn].size(); i++) {
			printf("%lf%c", ps[pn][i], (i == ps[pn].size() - 1) ? '\n' : ' ');
		}
	}
}B;
double x[BM::MAXN];
int main() {
	for (int n; ~scanf("%d", &n); ) {
		for (int i = 1; i <= n; i++) {
			scanf("%lf", &x[i]);
		}
		B.Solve(x, n), B.print();
	}
}
