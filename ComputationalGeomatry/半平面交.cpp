//对于给出点的顺时针和逆时针顺序不同,只需要加个 reverse 函数将点的顺序颠倒
int sgn(double x)
{ //符号函数
	if(fabs(x) < eps) return 0;
	if(x < 0) return -1;
	else return 1;
}
struct Point
{ //点
	double x,y;
	Point(){}
	Point(double _x,double _y)
	{
		x = _x; y = _y;
	}
	Point operator -(const Point &b)const
	{
		return Point(x - b.x, y - b.y);
	}
	double operator ^(const Point &b)const
	{ //叉积
		return x*b.y - y*b.x;
	}
	double operator *(const Point &b)const
	{ //点积
		return x*b.x + y*b.y;
	}
};
struct Line
{ //向量
	Point s,e; //两点
	double k; //斜率
	Line(){}
	Line(Point _s,Point _e)
	{ //构造
		s = _s; e = _e;
		k = atan2(e.y - s.y,e.x - s.x);
	}
	Point operator &(const Line &b)const
	{ //求两直线交点
		Point res = s;
		double t = ((s - b.s)^(b.s - b.e))/((s - e)^(b.s - b.e));
		res.x += (e.x - s.x)*t;
		res.y += (e.y - s.y)*t;
		return res;
	}
};
Line Q[MAXN];
Point p[MAXN]; //记录最初给的点集
Line line[MAXN]; //由最初的点集生成直线的集合
Point pp[MAXN]; //记录半平面交的结果的点集
//半平面交，直线的左边代表有效区域
bool HPIcmp(Line a,Line b)
{ //直线排序函数
	if(fabs(a.k - b.k) > eps)return a.k < b.k; //斜率排序
	//斜率相同我也不知道怎么办
	return ((a.s - b.s)^(b.e - b.s)) < 0;
}
void HPI(Line line[], int n, Point res[], int &resn)
{ //line是半平面交的直线的集合 n是直线的条数 res是结果
//的点集 resn是点集里面点的个数
	int tot = n;
	sort(line,line+n,HPIcmp);
	tot = 1;
	for(int i = 1;i < n;i++)
		if(fabs(line[i].k - line[i-1].k) > eps) //去掉斜率重复的
			line[tot++] = line[i];
	int head = 0, tail = 1;
	Q[0] = line[0];
	Q[1] = line[1];
	resn = 0;
	for(int i = 2; i < tot; i++)
	{
		if(fabs((Q[tail].e-Q[tail].s)^(Q[tail-1].e-Q[tail-1].s)) < eps ||
		fabs((Q[head].e-Q[head].s)^(Q[head+1].e-Q[head+1].s)) < eps)
			return;
		while(head < tail && (((Q[tail]&Q[tail-1]) -
		line[i].s)^(line[i].e-line[i].s)) > eps)
			tail--;
		while(head < tail && (((Q[head]&Q[head+1]) -
		line[i].s)^(line[i].e-line[i].s)) > eps)
			head++;
		Q[++tail] = line[i];
	}
	while(head < tail && (((Q[tail]&Q[tail-1]) -
	Q[head].s)^(Q[head].e-Q[head].s)) > eps)
		tail--;
	while(head < tail && (((Q[head]&Q[head-1]) -
	Q[tail].s)^(Q[tail].e-Q[tail].e)) > eps)
		head++;
	if(tail <= head + 1) return;
	for(int i = head; i < tail; i++)
		res[resn++] = Q[i]&Q[i+1];
	if(head < tail - 1)
		res[resn++] = Q[head]&Q[tail];
}
double dist(Point a,Point b)
{ //两点间距离
    return sqrt((a-b)*(a-b));
}
void change(Point a,Point b,Point &c,Point &d,double p)
{ //将线段ab往左移动距离p,修改得到线段cd
    double len=dist(a,b);
    /*三角形相似推出下面公式*/
    double dx=(a.y-b.y)*p/len;
    double dy=(b.x-a.x)*p/len;
    c.x=a.x+dx; c.y=a.y+dy;
    d.x=b.x+dx; d.y=b.y+dy;
}
double BSearch()
{ //二分搜索
	double l=0,r=100000;
    double ans=0;
    while(r-l>=eps)
    {
        double mid=(l+r)/2;
        for(int i=0;i < n;i++)
        {
            Point t1,t2;
            change(p[i],p[(i+1)%n],t1,t2,mid);
            line[i]=Line(t1,t2);
        }
        int resn;
        HPI(line,n,pp,resn);
        //等于0说明移多了
        if(resn==0) r=mid-eps;
        else l=mid+eps;
    }
    return l;
}
