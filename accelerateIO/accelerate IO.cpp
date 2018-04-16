#include<cstdio>
typedef long long ll;
inline ll in()
{
	char c=getchar();
	while(!isdigit(c))c=getchar();
	int x=0;
	while(isdigit(c))
	{
		x=x*10+c-'0';
		c=getchar();
	}
	return x;
}
void out(ll x)
{
	if(x>9)out(x/10);
	putchar(x%10+'0');
}
#define in in()
