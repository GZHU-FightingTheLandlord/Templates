/*
  Author: SemonChan
  Problem: BZOJ - 3884
  Problem Source: https://cn.vjudge.net/problem/HYSBZ-3884
                  https://www.lydsy.com/JudgeOnline/problem.php?id=3884
  Solution Source: https://cn.vjudge.net/solution/13601290
*/
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ll;

ll eular(ll n)
{
	ll res = n, a = n;
	for (ll i = 2;i * i <= a; i++)
	{
		if (a % i == 0)
		{
			res = res / i * (i - 1);
			while (a % i == 0)
				a /= i;
		}
	}
	if (a > 1)
		res = res / a * (a - 1);
	return res;
}


ll fpow(ll a, ll k, ll m)
{
	ll base = 1;
	if (a >= m)
		a = a % m + m;
	while (k)
	{
		if (k & 1)
		{
			base = base * a;
			if (base >= m)
				base = base % m + m;
		}
		a = a * a;
		if (a >= m)
			a = a % m + m;
		k >>= 1;
	}
	return base;
}


ll solve(ll m)
{
	if (m == 1)
		return 0;
	ll eu = eular(m);
	return fpow(2, solve(eu) + eu, m);
}

int main()
{
	ll n, m;
	scanf("%lld", &n);
	while (n--)
	{
		scanf("%lld", &m);
		printf("%lld\n", solve(m) % m);
	}
	return 0;
}
