/*
 * Source : http://172.22.27.1/problem?pid=1055
 * Author : ConanYu
 */

#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
ll a[500005];
int main()
{
	int T,n;
	ll m;
	scanf("%d",&T);
	while(T--)
	{
		memset(a,0x3f,sizeof(a));
		scanf("%d",&n);
		int cnt=1;
		scanf("%lld",&m);
		a[0]=m;
		for(int i=1;i<n;i++)
		{
			scanf("%lld",&m);
			if(m>=a[cnt-1])a[cnt++]=m;
			else {
				int pos=upper_bound(a,a+cnt,m)-a;
				a[pos]=m;
			}
		}
		printf("%d\n",n-cnt);
	}
}
