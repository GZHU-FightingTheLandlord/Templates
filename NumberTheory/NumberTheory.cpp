namespace NumberTheory
{
	typedef long long ll;
	ll ___inv,___t;
	void ex_gcd_fun(ll a, ll b,ll & ___inv,ll & ___t)
	{
		if(!b)
		{
			___inv=1;
			___t=0;
		}
		else
		{
			ex_gcd_fun(b,a%b,___t,___inv);
			___t-=___inv*(a/b);
		}
	}
	//above are private

	ll ex_gcd(ll a,ll b)
	{
		ex_gcd_fun(a,b,___inv,___t);
		return ___inv;
	}
	//ex_gcd(a,b) 求的是a在模b下的逆元
	//注意： 逆元有可能为负数

	ll eular(ll x)
	{
		//if(x<=0)return std::runtime_error("Eular error!");
		ll ans=x;
		for(int i=2;i*i<=x;i++)
		{
			if(x%i==0)
			{
				ans-=ans/i;
				while(x%i==0)x/=i;
			}
		}
		if(x!=1)ans-=ans/x;
		return ans;
	}

	const ll _MOD=1e9+7;
	//there change the MOD.
	ll Pow(ll a, ll b)
	{
		ll ans=1;
		a%=_MOD;
		while(b)
		{
			if(b&1)ans=(ans*a)%_MOD;
			a=(a*a)%_MOD;
			b>>=1;
		}
		return ans;
	}

	ll gcd(ll a,ll b)
	{
		return b==0?a:gcd(b,a%b);
	}

	ll lcm(ll a,ll b)
	{
		return a/gcd(a,b)*b;
	}

	bool isprime(ll x)
	{
		for(ll i=2;i*i<=x;i++)
		{
			if(x%i==0)return false;
		}
		return true;
	}

	void prime_table(bool *a,int len)
	{
		memset(a,0,len);
		a[0]=a[1]=false;
		for(int i=2;i<len;i++)if(!a[i])for(int j=i+i;j<len;j+=i)a[j]=true;
	}
	
	const LL MOD = 1e9 + 7;
	LL fun(LL a, LL b)
	{
		LL sum = 0;
		LL k = 1;
		while(b)
		{
			if(b&1)
				sum = (sum + a * k) % MOD;
			k = (k * 2) % MOD;
			b >>= 1;
		}
		return sum;
	}

	//it will be updated
}
