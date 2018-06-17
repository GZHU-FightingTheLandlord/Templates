#include<cstdio>
#include<cctype>
struct InputStream
{
	const static int BUF_SIZE = 100000;
#define llong long long
	inline char nc()
	{
		static char buf[BUF_SIZE], *p1 = buf + BUF_SIZE, *pend = buf + BUF_SIZE;
		if(p1 == pend)
		{
			p1 = buf;
			pend = buf + fread(buf, 1, BUF_SIZE, stdin);
		}
		return (*p1++);
	}
	inline llong in()
	{
		char c = nc();
		while(!isdigit(c))
			c = nc();
		llong x = 0;
		while(isdigit(c))
		{
			x = x * 10 + c - '0';
			c = nc();
		}
		return x;
	}
	template<typename T>
	inline void g(T & x)
	{
		x = in();
	}
	inline llong g()
	{
		return in();
	}
#undef llong
} in;
#define rd in.g()

//int a = rd;
//in.g(a);
