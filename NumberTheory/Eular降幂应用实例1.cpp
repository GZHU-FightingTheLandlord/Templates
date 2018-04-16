/*
    Author: SemonChan
    Problem: eoj-3394 瓜皮快速幂
    Problem Source: https://acm.ecnu.edu.cn/problem/3394/
*/
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

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

ll mul(ll a, ll b, ll m)
{
	return (a * b - ll(a / (long double)m * b + 1e-3) * m + m) % m;
}

ll fpow(ll a, ll k, ll m)
{
	ll base = 1;
	while (k)
	{
		if (k & 1)
			base = mul(base, a, m);
		a = mul(a, a, m);
		k >>= 1;
	}
	return base;
}

ll solve(ll n, ll m)
{
	if (n == 0 || m == 1)
		return 2 % m;
	if (n == 1)
		return 4 % m;
	if (n == 2)
		return 16 % m;
	if (n == 3)
		return 65536 % m;
   	ll eu = eular(m);
	return fpow(2, solve(n - 1, eu) + eu, m);
}

int main()
{
	ll n, m;
	scanf("%lld%lld", &n, &m);
	printf("%lld\n", solve(n, m));
	return 0;
}
