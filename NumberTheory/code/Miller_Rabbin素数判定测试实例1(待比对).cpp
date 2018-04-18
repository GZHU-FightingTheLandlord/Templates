/*
 * check()	参数: num
 *
 * Time Complexity: O(logN)
 *
 * 错误概率: 2^(-k)
 *
 * 直接指定方式: n<1e6, k=2, a=2/a=3; n为32位无符号整数, k=3,a=2/a=7/a=61
 */

struct Miller_Rabbin
{
	LL MyRand(LL e)
	{
		return rand() % e + 1;
	}

	LL MultMod(LL a, LL b, LL c)
	{
		a %= c, b %= c;
		LL ans = 0;
		while(b)
		{
			if(b & 1)
				ans += a, ans %= c;
			a <<= 1;
			if(a >= c)
				a %= c;
			b >>= 1;
		}
		return ans;
	}

	LL gt(LL a, LL u, LL num)
	{
		LL cur = 1, nxt = a;
		while(u)
		{
			if((u & 1) > 0)
				cur = MultMod(cur, nxt, num);
			else
				nxt = MultMod(nxt, nxt, num);
			u = u>>1;
		}
		return cur % num;
	}

	bool check(LL num)
	{
		if(num == 2)
			return true;
		if(num < 2 || num % 2 == 0)
			return false;
		LL u = num - 1;
		int times = 20;	//checking times
		while(u % 2 == 0)
			u /= 2;
		for(int i = 0; i < times; i++)
		{
			LL a = MyRand(num - 1);
			LL x = gt(a, u, num);
			LL y = x;
			LL tu = u;
			while(tu < num)
			{
				y = MultMod(x, x, num);
				if(y == 1 && x != 1 && x != num - 1)
					return false;
				x = y;
				tu *= 2;
			}
			if(x != 1)
				return false;
		}
		return true;
	}
};
