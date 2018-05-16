/*
	Author: Seee
	Problem Source: POJ-2299
	Code Source: https://cn.vjudge.net/solution/13941687
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <bitset>
#include <math.h>
#include <cmath>
#include <time.h>
#include <string.h>
#include <vector>
#include <set>
#include <deque>
#include <stack>
#include <time.h>
#include <map>
#include <queue>
#include <functional>
#include <cctype>
#include <iomanip>

#pragma warning(disable:4996)
#pragma comment(linker, "/STACK:336777216")
using namespace std;

#define pb(a)        push_back(a)
#define mp(a, b)     make_pair(a, b)
#define all(a)       a.begin(),a.end()
#define szz(a)       (int)a.size()
#define endl         '\n'

typedef long long ll;
typedef double db;
typedef long double ld;
typedef unsigned long long ull;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef pair<ll, int> pli;
typedef pair<db, db> pdd;

const ll MOD = 1e9 + 7;
const int inf = 0x3f3f3f3f;
const int INF = 0x7fffffff;
const ll llINF = 0x3f3f3f3f3f3f3f3f;
const ll LLINF = 0x7fffffffffffffff;
const db PI = acos(-1.0);
const db eps = 1e-8;

inline ll g() { char c; ll p = 1; while (!isdigit(c = getchar())) if (c == '-') p = -p; ll x = c - '0'; while (isdigit(c = getchar())) x = x * 10 + c - '0'; return p * x; }

class FenwickTree {
public:
	void update(int i, int delta)
	{
		while (i <= n) {
			sum_[i] += delta;
			i += lowbit(i);
		}
	}

	int query(int l, int r)
	{
		return getres(r) - getres(l - 1);
	}

	FenwickTree(int _n = 0) : sum_(_n + 1, 0), n(_n) {}
private:
	int n;
	vector<int> sum_;

	int getres(int i)
	{
		int res = 0;
		while (i > 0) {
			res += sum_[i];
			i -= lowbit(i);
		}
		return res;
	}

	inline int lowbit(int x)
	{
		return x & (-x);
	}
};

int v[600005];
int tmp[600005];
int main()
{
	int n;
	while (scanf("%d", &n) == 1 && n) {
		for (int i = 0; i < n; i++) scanf("%d", &v[i]), tmp[i] = v[i];

		/* 离散	*/
		sort(tmp, tmp + n);
		int m = unique(tmp, tmp + n) - tmp;
		for (int i = 0; i < n; i++) {
			int pos = lower_bound(tmp, tmp + m, v[i]) - tmp;
			v[i] = pos + 1;
		}


		FenwickTree f(n);
		ll ans = 0;
		for (int i = n - 1; i >= 0; i--) {
			ans += f.query(1, v[i]);
			f.update(v[i] + 1, 1);
		}
		printf("%lld\n", ans);
	}
}
