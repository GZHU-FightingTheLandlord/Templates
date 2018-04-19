/*
 * setvalue() 参数: n设置种类数
 * init() 无参数,进行预处理
 * solve() 参数: n数值
 *
 * Time Complexity: O(n^3)
 *
 */


struct GeneratingFunction
{
	LL c1[105], c2[105];
	int N;

	void setvalue(int n)
	{
		N = n;
	}

	void init(void)
	{
		c1[0] = 1;
		for(int i = 1; i <= N; ++i)
		{
			for(int j = 0; j < 100; j += i)
				for(int k = 0; j + k < 100; ++k)
					c2[j + k] += c1[k];
			for(int j = 0; j < 100; ++j)
				c1[j] = c2[j], c2[j] = 0;
		}
	}

	LL solve(int n)
	{
		return c1[n];
	}
};
