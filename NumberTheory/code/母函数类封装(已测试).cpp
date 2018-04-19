/*
 *
 * init(): 参数: 所取最大值, 调用进行预处理打表
 * solve(): 返回n组合的最大值
 *
 * Time Complexity: O(n^3)
 *
 * 转换公式: (1+x^2+x^3+x^4+...+x^n)*(1+x^2+x^4+x^6+x^8+...)*(1+x^3+x^6+x^9+...) --- 共有n个表达式相乘
 *
 * 讲解: https://blog.csdn.net/HowardEmily/article/details/75041523
 *
 * 应用实例: hdu1028 / hdu1398
 *
 */

struct GeneratingFunction
{
	long long pre[1005], now[1005];

	void init(int N)
	{
		memset(pre, 0, sizeof(pre));
		memset(now, 0, sizeof(now));
		for(int i = 0; i <= N; ++i)
			pre[i] = 1;

		for(int i = 2; i <= N; ++i)
		{
			for(int j = 0; j <= N; ++j)
			{
				for(int k = 0; j + k <= N; k += i)
				{
					now[j + k] += pre[j];
				}
			}
			for(int j = 0; j <= N; j++)
				pre[j] = now[j], now[j] = 0;
		}
	}

	long long solve(int n)
	{
		return pre[n];
	}
};
