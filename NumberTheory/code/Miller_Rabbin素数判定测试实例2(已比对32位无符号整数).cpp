/*
	IsPrime(): 参数: n -- 需要判断的素数
	Time Complexity: O(1)
 	错误概率: 2^(-k)
	直接指定方式: n<1e6, k=2, a=2/a=3; n为32位无符号整数, k=3,a=2/a=7/a=61 
*/
struct Miller_Rabbin
{
	LL fpow(LL a, LL b, LL r)
	{
	    LL ans=1,buff=a;
	    while(b)
	    {
	        if(b&1)
	            ans = (ans * buff) % r;
	        buff = (buff * buff) % r;
	        b>>=1;
	    }
	    return ans;
	}
	bool Miller(LL n, LL a)
	{
	    LL r = 0, s = n - 1, j;
	    if(!(n % a))
	        return false;
	    while(!(s&1))
	        s >>= 1, r++;
	    LL k = fpow(a, s, n);
	    if(k == 1)
	        return true;
	    for(j = 0; j < r; j++, k = k * k % n)
	        if(k == n - 1)
	            return true;
	    return false;
	}
	bool IsPrime(int n)
	{
		if(n < 2)
			return false;
	    int tab[] = {2, 3, 5, 7};
	    for(int i = 0; i < 4; i++)
	    {
	        if(n == tab[i])
	            return true;
	        if(!Miller(n, tab[i]))
	            return false;
	    }
	    return true;
	}
};
