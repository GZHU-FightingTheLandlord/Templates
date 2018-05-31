#include <bits/stdc++.h>
#pragma comment(linker,"/STACK:102400000,102400000")
using namespace std;
const int MAXN=10005;
class berlekamp_massey
{
private:
	int n;
	vector<double>ps[MAXN];
	int pn, fail[MAXN];
	double delta[MAXN];
public:
	void solve(double *x, int n)
	{
		pn=0;
		memset(fail, 0, sizeof(fail));
		memset(delta, 0, sizeof(delta));
		ps[0].clear();
		for(int i = 1; i < n + 1; i++)
		{
			double dt=-x[i];
			for(int j = 0; j < ps[pn].size(); j++)
			{
				dt += x[i - j - 1] * ps[pn][j];
			}
			delta[i] = dt;
			if(fabs(dt) <= 1e-8)continue;
			fail[pn] = i;
			if(!pn)
			{
				ps[++pn].resize(1);
				continue;
			}
			vector<double>& ls = ps[pn - 1];
			double k = -(dt / delta[fail[pn - 1]]);
			vector<double> cur;
			cur.resize(i - fail[pn - 1] - 1);
			cur.push_back(-k);
			for(int j = 0; j < ls.size(); j++)
			{
				cur.push_back(ls[j] * k);
			}
			if(cur.size() < ps[pn].size())
			{
				cur.resize(ps[pn].size());
			}
			for(int j = 0; j < ps[pn].size(); j++)
			{
				cur[j] += ps[pn][j];
			}
			ps[++pn] = cur;
		}
	}
	void print()
	{
		for(int g = 0; g < ps[pn].size(); g++)
		{
			printf("%lf ", ps[pn][g]);
		}
		printf("\n");
	}
}solver;

double x[MAXN];

int main()
{
	int n;
	while(~scanf("%d", &n))
	{
		for(int i = 1; i < n + 1; i++)
		{
			scanf("%lf", x + i);
		}
		solver.solve(x, n);
		solver.print();
	}
}
