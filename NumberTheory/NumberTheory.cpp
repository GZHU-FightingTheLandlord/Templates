#include<cstdio>
typedef long long ll;
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
	//ex_gcd(a,b) �����a��ģb�µ���Ԫ
	//ע�⣺ ��Ԫ�п���Ϊ����

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

	//it will be updated
}

int main()
{
	using namespace NumberTheory;
	printf("%lld\n",Pow(15151,641654));
}
