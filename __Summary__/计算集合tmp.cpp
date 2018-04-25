#include <bits/stdc++.h>

using namespace std;

typedef struct Point Point;
typedef struct Line Line;
typedef struct Polygon Polygon;
typedef struct Polygon_convex Polygon_convex;
typedef struct Halfline Halfline;
typedef complex<double> point;
typedef pair<point, point> halfplane;

const double eps = 1e-8;
const double pi = acos(-1.0);
const int maxN = 1e3 + 5;

bool comp_less(const Point &a, const Point &b);
Polygon_convex convex_hull(vector<Point> a);
bool containOn(const Polygon_convex &a, const Point &b);	//判断点是否在凸包内
int containOlogn(const Polygon_convex &a, const Point &b);	//求出凸包
double convex_diameter(Polygon_convex &a, int &First, int &Second);	//求凸多边形的直径
inline double sqr(double x);	//计算一个数的平方
int cmp(double x);	//计算几何误差修正
int gcd(int a, int b);
double det(const Point &a, const Point &b);	//计算两个向量的叉积
double dot(const Point &a, const Point &b);	//计算两个向量的点积
double dist(const Point &a, const Point &b);	//计算两个点的距离
Point rotate_point(const Point &p, double A);	//op绕远点逆时针旋转A(弧度)
Line point_make_line(const Point a, const Point b);	//用两个点a,b生成的一个线段或者直线
double dis_point_segment(const Point p, const Point s, const Point t);	//求p点到线段st的距离
void PointProjLine(const Point p, const Point s, const Point t, Point &cp);	//求p点到线段st的垂足,保存在cp中
bool PointOnSegment(Point p, Point s, Point t);	//判断p点是否在线段st上(包括端点)
bool parallel(Line a, Line b);	//判断a和b是否平行
bool line_make_point(Line a, Line b, Point &res);	//判断a和b是否相交,如果相交则返回true且交点保存在res中
Line move_d(Line a, const double &len); //将直线a沿法向量方向平移距离len得到的直线

struct Point
{
	double x, y;
	Point(){}
	Point(double a, double b):x(a), y(b){}
	void input()
	{
		scanf("%lf%lf", &x, &y);
	}
	friend Point operator +(const Point &a, const Point &b)
	{
		return Point(a.x + b.x, a.y + b.y);
	}
	friend Point operator -(const Point &a, const Point &b)
	{
		return Point(a.x - b.x, a.y - b.y);
	}
	friend bool operator ==(const Point &a, const Point &b)
	{
		return cmp(a.x - b.x) == 0 && cmp(a.y - b.y) == 0;
	}
	friend Point operator *(const Point &a, const double &b)
	{
		return Point(a.x * b, a.y * b);
	}
	friend Point operator *(const double &a, const Point &b)
	{
		return Point(a * b.x, a * b.y);
	}
	friend Point operator /(const Point &a, const double &b)
	{
		return Point(a.x / b, a.y / b);
	}
	double norm()
	{
		return sqrt(sqr(x) + sqr(y));
	}
};

struct Line
{
	Point a, b;
	Line() {}
	Line(Point x, Point y):a(x), b(y) {}
};

struct Polygon
{
	int n;
	Point a[maxN];
	Polygon() {}
	double perimeter()	//计算多边形周长
	{
		double sum = 0;
		a[n] = a[0];
		for(int i = 0; i < n; i++)
			sum += (a[i + 1] - a[i]).norm();
		return sum;
	}
	double area()	//计算多边形面积
	{
		double sum = 0;
		a[n] = a[0];
		for(int i = 0; i < n; i++)
			sum += det(a[i + 1], a[i]);
		return sum / 2;
	}
	int Point_In(Point t)	//判断点是否在多边形内部
	{
		int num = 0;
		a[n] = a[0];
		for(int i = 0; i < n; i++)
		{
			if(PointOnSegment(t, a[i], a[i + 1]))
				return 2;
			int k = cmp(det(a[i + 1] - a[i], t - a[i]));
			int d1 = cmp(a[i].y - t.y);
			int d2 = cmp(a[i + 1].y - t.y);
			if(k > 0 && d1 <= 0 && d2 > 0)
				num++;
			if(k < 0 && d2 <= 0 && d1 > 0)
				num--;
		}

	}
	Point MassCenter()	//求多边形的重心
	{
		Point ans = Point(0, 0);
		if(cmp(area()) == 0)
			return ans;
		a[n] = a[0];
		for(int i = 0; i < n; i++)
			ans = ans + (a[i] + a[i + 1]) * det(a[i + 1], a[i]);
		return ans / area() / 6.;
	}
	int Border_Int_Point_Num()	//求多边形边界上的格点个数
	{
		int num = 0;
		a[n] = a[0];
		for(int i = 0; i < n; i++)
			num += gcd(fabs(int(a[i + 1].x - a[i].x)), fabs(int(a[i + 1].y - a[i].y)));
		return num;
	}
	int Inside_Int_Point_Num()	//求多边形内的格点个数
	{
		return int(area()) + 1 - Border_Int_Point_Num() / 2;
	}
};

struct Polygon_convex
{
	vector<Point> P;
	Polygon_convex(int Size = 0)
	{
		P.resize(Size);
	}
};

struct Halfline
{
	double a, b, c;
	Halfline(Point p, Point q)
	{
		a = p.y - q.y;
		b = q.x - p.x;
		c = det(p, q);
	}
	Halfline(double aa, double bb, double cc)
	{
		a = aa, b = bb, c = cc;
	}
	double calc(Halfline &L, Point &a)
	{
		return a.x * L.a + a.y * L.b + L.c;
	}
	Point Intersect(Point &a, Point &b, Halfline &L)
	{
		Point res;
		double t1 = calc(L, a), t2 = calc(L, b);
		res.x = (t2 * a.x - t1 * b.x) / (t2 - t1);
		res.y = (t2 * a.y - t1 * b.y) / (t2 - t1);
		return res;
	}
	Polygon_convex cut(Polygon_convex &a, Halfline &L)
	{
		int n = a.P.size();
		Polygon_convex res;
		for(int i = 0; i < n; ++i)
		{
			if(calc(L, a.P[i]) < -eps)
				res.P.push_back(a.P[i]);
			else
			{
				int j = i - 1;
				if(j < 0)
					j = n - 1;
				if(calc(L, a.P[j]) < -eps)
					res.P.push_back(Intersect(a.P[j], a.P[i], L));
				j = i + 1;
				if(j == n)
					j = 0;
				if(calc(L, a.P[j]) < -eps)
					res.P.push_back(Intersect(a.P[i], a.P[j], L));
			}
		}
		return res;
	}
};

bool comp_less(const Point &a, const Point &b)
{
	return cmp(a.x - b.x) < 0 || cmp(a.x - b.x) == 0 && cmp(a.y - b.y) < 0;
}

Polygon_convex convex_hull(vector<Point> a)
{
	Polygon_convex res(2 * a.size() + 5);
	sort(a.begin(), a.end(), comp_less);
	a.erase(unique(a.begin(), a.end()), a.end());
	int m = 0;
	for(int i = 0; i < (int)a.size(); ++i)
	{
		while(m > 1 && cmp(det(res.P[m - 1] - res.P[m - 2], a[i] - res.P[m - 2])) <= 0)
			--m;
	}
	int k = m;
	for(int i = int(a.size()) - 2; i >= 0; --i)
	{
		while(m > k && cmp(det(res.P[m - 1] - res.P[m - 2], a[i] - res.P[m - 2])) <= 0)
			--m;
		res.P[m++] = a[i];
	}
	res.P.resize(m);
	if(a.size() > 1)
		res.P.resize(m - 1);
	return res;
}

bool containOn(const Polygon_convex &a, const Point &b)
{
	int n = a.P.size();
	#define next(i) ((i + 1) % n)
	int sign = 0;
	for(int i = 0; i < n; ++i)
	{
		int x = cmp(det(a.P[i] - b, a.P[next(i)] - b));
		if(x)
		{
			if(sign)
				if(sign != x)
					return false;
			else
				sign = x;
		}
	}
	return true;
}

int containOlogn(const Polygon_convex &a, const Point &b)
{
	int n = a.P.size();
	//找一个凸包内部的点g
	Point g = (a.P[0] + a.P[n / 3] + a.P[2 * n / 3]) / 3.;
	int l = 0, r = n;
	//二分凸包g - a.P[a] - a.P[b]
	while(l + 1 < r)
	{
		int mid = (l + r) / 2;
		if(cmp(det(a.P[l] - g, a.P[mid] - g)) > 0)
		{
			if(cmp(det(a.P[l] - g, b - g)) >= 0 && cmp(det(a.P[mid] - g, b - g)) < 0)
				r = mid;
			else
				l = mid;
		}
		else
		{
			if(cmp(det(a.P[l] - g, b - g)) < 0 && cmp(det(a.P[mid] - g, b - g)) >= 0)
				l = mid;
			else
				r = mid;
		}
	}
	r %= n;
	int z = cmp(det(a.P[r] - b, a.P[l] - b)) - 1;
	if(z == -2)
		return 1;
	return z;

}

double convex_diameter(Polygon_convex &a, int &First, int &Second)
{
	vector<Point> &p = a.P;
	int n = p.size();
	double maxd = 0.;
	if(n == 1)
	{
		First = Second = 0;
		return maxd;
	}
	#define next(i) ((i + 1) % n)
	for(int i = 0, j = 1; i < n; ++i)
	{
		while(cmp(det(p[next(i)] - p[i], p[j] - p[i]) - det(p[next(i)] - p[i], p[next(j)] - p[i])) < 0)
			j = next(j);
		double d = dist(p[i], p[j]);
		if(d > maxd)
			maxd = d, First =i, Second = j;
		d = dist(p[next(i)], p[next(j)]);
		if(d > maxd)
			maxd = d, First = i, Second = j;
	}
	return maxd;
}

inline double sqr(double x)
{
	return x * x;
}

int gcd(int a, int b)
{
	return b == 0 ? a : gcd(b, a % b);
}

int cmp(double x)
{
	if(fabs(x) < eps)
		return 0;
	if(x > 0)
		return 1;
	return -1;
}

double det(const Point &a, const Point &b)
{
	return a.x * b.y - a.y * b.x;
}

double dot(const Point &a, const Point &b)
{
	return a.x * b.x + a.y * b.y;
}

double dist(const Point &a, const Point &b)
{
	return (a - b).norm();
}

Point rotate_point(const Point &p, double A)
{
	double tx = p.x, ty = p.y;
	return Point(tx * cos(A) - ty * sin(A), tx * sin(A) + ty * cos(A));
}

Line point_make_line(const Point a, const Point b)
{
	return Line(a, b);
}

double dis_point_segment(const Point p, const Point s, const Point t)
{
	if(cmp(dot(p - s, t - s)) < 0)
		return (p - s).norm();
	if(cmp(dot(p - t, s - t)) < 0)
		return (p - t).norm();
	return fabs(det(s - p, t - p) / dist(s, t));
}

void PointProjLine(const Point p, const Point s, const Point t, Point &cp)
{
	double r = dot((t - s), (p - s)) / dot(t - s, t - s);
	cp = s + r * (t - s);
}

bool PointOnSegment(Point p, Point s, Point t)
{
	return cmp(det(p - s, t - s)) == 0 && cmp (dot(p - s, p - t)) <= 0;
}

bool parallel(Line a, Line b)
{
	return !cmp(det(a.a - a.b, b.a - b.b));
}

bool line_make_point(Line a, Line b, Point &res)
{
	if(parallel(a, b))
		return false;
	double s1 = det(a.a - b.a, b.b - b.a);
	double s2 = det(a.b - b.a, b.b - b.a);
	res = (s1 * a.b - s2 * a.a) / (s1 - s2);
	return true;
}

Line move_d(Line a, const double &len)
{
	Point d = a.b - a.a;
	d = d / d.norm();
	d = rotate_point(d, pi / 2);
	return Line(a.a + d * len, a.b + d * len);
}

int main()
{
	return 0;
}
