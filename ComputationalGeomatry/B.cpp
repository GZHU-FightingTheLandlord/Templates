#include<cstdio>
#include<cmath>
#include<algorithm>
#include<cstdlib>
#include<ctime>
using namespace std;
struct point
{
	double x;
	double y;
	point(double _x=0.0,double _y=0.0):x(_x),y(_y){}
}p[1005],O;
double R;
const double eps=1e-7;
double sqr(const double & x)
{
	return x*x;
}
double dis(const point & x,const point & y)
{
	return sqrt(sqr(x.x-y.x)+sqr(x.y-y.y));
}
bool isin(const point & x)
{
	if(dis(O,x)<=R+eps)return true;
	return false;
}
point midd(const point & x,const point & y)
{
	double a=(x.x+y.x)/2.0;
	double b=(x.y+y.y)/2.0;
	return point(a,b);
}
struct line
{
	point a,b;
	line(){}
	line(point _a,point _b):a(_a),b(_b){}
};
point cross(const line & u,const line & v)
{
	point ret=u.a;
	double t1=(u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x);
	double t2=(u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x);
	double t=t1;
	t/=t2;
	ret.x+=(u.b.x-u.a.x)*t;
	ret.y+=(u.b.y-u.a.y)*t;
	return ret;
}
point cir(const point & a,const point & b,const point & c)
{
	line u,v;
	u.a.x=(a.x+b.x)/2;
	u.a.y=(a.y+b.y)/2;
	u.b.x=u.a.x-a.y+b.y;
	u.b.y=u.a.y+a.x-b.x;
	v.a.x=(a.x+c.x)/2;
	v.a.y=(a.y+c.y)/2;
	v.b.x=v.a.x-a.y+c.y;
	v.b.y=v.a.y+a.x-c.x;
	return cross(u,v);
}
int main()
{
	srand(time(0));
	int X,Y,n;
	while(~scanf("%d %d %d",&X,&Y,&n))
	{
		for(int i=0;i<n;i++)scanf("%lf %lf",&p[i].x,&p[i].y);
		random_shuffle(p,p+n);
		O.x=O.y=R=0.0;
        for(int i=0;i<n;i++)
		{
			if(!isin(p[i]))
			{
				O=p[i];R=0.0;
				for(int j=0;j<i;j++)
				{
					if(!isin(p[j]))
					{
						O=midd(p[i],p[j]);
						R=dis(O,p[j]);
						for(int k=0;k<j;k++)
						{
							if(!isin(p[k]))
							{
								O=cir(p[i],p[j],p[k]);
								R=dis(O,p[i]);
							}
						}
					}
				}
			}
		}
		printf("(%.1lf,%.1lf).\n",O.x,O.y);
		printf("%.1lf\n",R);
	}
}
