#include<cstdlib>
#include<cmath>
#include<ctime>
#include<iostream>
#include<iomanip>
#include<algorithm>
using namespace std;

struct point
{
	double x;
	double y;
	point(){}
	point(double _x,double _y):x(_x),y(_y){}
}e[105];

int n;

double sqr(const double & x)
{
	return x*x;
}

double dis(const point & x,const point & y)
{
	return sqrt(sqr(x.x-y.x)+sqr(x.y-y.y));
}

double cal(const point & x)
{
	double ans=0.0;
	for(int i=0;i<n;i++)
	{
		ans+=dis(x,e[i]);
	}
	return ans;
}

point p;

point sca(const double & t)
{
	double x=p.x+(rand()&1?1:-1)*rand()*t/RAND_MAX;
	double y=p.y+(rand()&1?1:-1)*rand()*t/RAND_MAX;
	return point(x,y);
}

const int TIMES=66;
const double ACCEPT_RATE=1e5;
const double M_TEMPERATRUE=0.93;

double SA()
{
	double mit=cal(p);
    int tms=0; // times
	for(double t=5000;t>1&&tms<TIMES;t*=M_TEMPERATRUE,tms++)
	{
		const point tp=sca(t);
		const double tmp=cal(tp);
		if(tmp<mit||rand()*ACCEPT_RATE<t*RAND_MAX)
		{
            p=tp;
            mit=tmp;
            tms=0;
		}
	}
	return mit;
}

int main()
{
	srand(unsigned(&n)*CLOCKS_PER_SEC/RAND_MAX);

	cin>>n;
	for(int i=0;i<n;i++)cin>>e[i].x>>e[i].y;

	p=e[rand()%n];

    cout<<(long long)(SA()+0.5)<<endl;
}
