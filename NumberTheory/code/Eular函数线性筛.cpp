/*
 * solve(n): 打表1~n内的欧拉函数和素数, 输出到euler[]和prime[]内
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
		isprime[0] = isprime[1] = false;
		euler[1] = 1;
	}

	void solve(int n)
	{
		init();

		for(int i = 2; i <= n; ++i)
		{
			if(isprime[i])
				prime[++tot] = i, euler[i] = i - 1;
			for(int j = 1; j <= tot && i * prime[j] <= n; j++)
			{
				int tmp = i * prime[j];
				isprime[tmp] = false;
				if(i % prime[j] == 0)
				{
					euler[tmp] = euler[i] * prime[j];
					break;
				}
				else
					euler[tmp] = euler[tmp] * (prime[j] - 1);
			}
		}
	}
};
