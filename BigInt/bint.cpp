/*
  Author : Caproner
*/

struct bint
{
	int a[1005];
	int n;
	
	bint()
	{
		n=1;
		memset(a,0,sizeof(a));
	}
	
	void init(char *s,int len)
	{
		n=len;
		for(int i=0;i<len;i++)
			a[i]=(s[len-i-1]-'0');
	}
	
	void mul_10()
	{
		for(int i=n;i>0;i--)
			a[i]=a[i-1];
		a[0]=0;
		n++;
	}
	
	bint operator +(const bint &b)const
	{
		bint c;
		c.n=0;
		int d=0;
		for(int i=0;d||i<n||i<b.n;i++)
		{
			c.a[i]=a[i]+b.a[i]+d;
			d=c.a[i]/10;
			c.a[i]%=10;
			c.n++;
		}
		return c;
	}
	
	bint operator -(const bint &b)const
	{
		bint c;
		c.n=0;
		int d=0;
		for(int i=0;d||i<n||i<b.n;i++)
		{
			c.a[i]=a[i]-b.a[i]-d;
			d=0;
			if(c.a[i]<0)
			{
				d=1;
				c.a[i]+=10;
			}
			c.n++;
		}
		while(c.n>1&&c.a[c.n-1]==0)
			c.n--;
		return c;
	}
	
	bint operator *(const bint &b)const
	{
		bint c;
		bint mul[10];
		for(int i=1;i<10;i++)
			mul[i]=mul[i-1]+(*this);
		for(int i=0;i<b.n;i++)
		{
			c=c+mul[b.a[i]];
			for(int j=0;j<10;j++)
				mul[j].mul_10();
		}
		return c;
	}
	
	void print()
	{
		for(int i=n-1;i>=0;i--)
			printf("%d",a[i]);
		printf("\n");
	}
	
}; 
