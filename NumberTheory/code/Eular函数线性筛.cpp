/*
 * solve(n): 打表1~n内的欧拉函数, 输出到ans[]内
 *
 * Time Complexity: O(n)
 *
 * P.S. 增加MAX值时可能爆栈, 手动加栈, 如下
 *
 * #pragma comment(linker, "/STACK:1024000000,1024000000")
 *
 */


struct GetEuler
{
	int prime[200000+5], euler[200000+5];
	bool isprime[200000+5];
	int tot;

	void init(void)
	{
		memset(isprime, true, sizeof(isprime));
		memset(prime, 0, sizeof(prime));
		memset(euler, 0, sizeof(euler));
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

	void solve(int n)
	{
		init();

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
