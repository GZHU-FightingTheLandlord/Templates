/*
 * 打表1~n内的欧拉函数和素数, 输出到Object.euler[], Object.prime[]内
 *
 * 参数: n
 * 
 * Time Complexity: O(n)
 *
 */

#include <bits/stdc++.h>

using namespace std;


struct GetEuler
{
	int *prime, *euler;
	bool *isprime;
	int tot;
	int len;

	GetEuler(int _len)
	{
		const int __len = _len * 2;
		len = _len * 2;
		prime = new int [__len];
		euler = new int [__len];
		isprime = new bool [__len];
		solve();
	}

	void init(void)
	{
		memset(isprime, true, sizeof(bool) * len);
		memset(prime, 0, sizeof(int) * len);
		memset(euler, 0, sizeof(int) * len);
		tot = 0;
	}

	void getprime(int n)
	{
		isprime[0] = isprime[1] = false;
		for(int i = 2; i <= n; i++)
		{
			if(isprime[i])
				prime[++tot] = i;
			for(int j = 1; j <= tot && i * prime[j] <= n; j++)
			{
				isprime[i * prime[j]] = false;
				if(i % prime[j] == 0)
					break;
			}
		}
	}

	void solve()
	{
		init();

		int n = len;

		getprime(2 * n);

		euler[1] = 1;
		for(int i = 2; i <= n; ++i)
		{
			if(isprime[i])
			{
				euler[i] = i - 1;
				continue;
			}
			for(int j = 1; j <= tot; ++j)
			{
				if(i % prime[j])
					continue;
				if(i / prime[j] % prime[j])
					euler[i] = euler[i / prime[j]] * (prime[j] - 1);
				else
					euler[i] = euler[i / prime[j]] * prime[j];
				break;
			}
		}
	}
};
