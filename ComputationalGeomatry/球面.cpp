#include <math.h>
const double pi=acos(-1);
//计算圆心角lat表示纬度,-90<=w<=90,lng表示经度
//返回两点所在大圆劣弧对应圆心角,0<=angle<=pi
double angle(double lng1,double lat1,double lng2,double lat2){
	double dlng=fabs(lng1-lng2)*pi/180;
	while (dlng>=pi+pi)
		dlng-=pi+pi;
	if (dlng>pi)
		dlng=pi+pi-dlng;
	lat1*=pi/180,lat2*=pi/180;
	return acos(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2));
}
//计算距离,r为球半径
double line_dist(double r,double lng1,double lat1,double lng2,double lat2){
	double dlng=fabs(lng1-lng2)*pi/180;
	while (dlng>=pi+pi)
		dlng-=pi+pi;
	if (dlng>pi)
		dlng=pi+pi-dlng;
	lat1*=pi/180,lat2*=pi/180;
	return r*sqrt(2-2*(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2)));
}
//计算球面距离,r为球半径
inline double sphere_dist(double r,double lng1,double lat1,double lng2,double lat2){
	return r*angle(lng1,lat1,lng2,lat2);
}
//球面反射
#include <cstdio>
#include <cmath>
const int size = 555;
const double eps = 1e-9;
struct point {double x, y, z;} centre = {0, 0, 0};
struct circle {point o; double r;} cir[size];
struct ray {point s, dir;} l;
int n;
int dcmp (double x){return x < -eps ? -1 : x > eps;}
double sqr (double x){return x*x;}
double dot (point a, point b){return a.x * b.x + a.y * b.y + a.z * b.z;}
double dis2 (point a, point b){return sqr(a.x-b.x) + sqr(a.y-b.y) + sqr(a.z-b.z);}
double disToLine2 (point a, ray l){/**** 点到直线L的距离的平方 **/
	point tmp;
	tmp.x =  l.dir.y * (a.z - l.s.z) - l.dir.z * (a.y - l.s.y);
	tmp.y = -l.dir.x * (a.z - l.s.z) + l.dir.z * (a.x - l.s.x);
	tmp.z =  l.dir.x * (a.y - l.s.y) - l.dir.y * (a.x - l.s.x);
	return dis2 (tmp, centre) / dis2 (l.dir, centre);
}
/**** 用向量法求交点  ***/
bool find (circle p, ray l, double &k, point &t)
{
	double h2 = disToLine2 (p.o, l);
//	printf ("h2 = %lf\n", h2);
	if (dcmp(p.r*p.r - h2) < 0) return false;
	point tmp;
	tmp.x = p.o.x - l.s.x;
	tmp.y = p.o.y - l.s.y;
	tmp.z = p.o.z - l.s.z;
	if (dcmp(dot(tmp, l.dir)) <= 0) return false;
	k = sqrt(dis2(p.o, l.s) - h2) - sqrt(p.r*p.r - h2);
	double k1 = k / sqrt(dis2(l.dir, centre));
	t.x = l.s.x + k1 * l.dir.x;
	t.y = l.s.y + k1 * l.dir.y;
	t.z = l.s.z + k1 * l.dir.z;
	return true;
}
/*计算新射线的起点和方向 */
void newRay (ray &l, ray l1, point inter)
{
	double k = - 2 * dot(l.dir, l1.dir);
	l.dir.x += l1.dir.x * k;
	l.dir.y += l1.dir.y * k;
	l.dir.z += l1.dir.z * k;
	l.s = inter;
}
/* 返回的是最先相交的球的编号,均不相交,返回-1 */
int update ()
{
	int sign = -1, i;
	double k = 1e100, tmp;
	point inter, t;
	for (i = 1; i <= n; i++){ //找到最先相交的球
		if (!find (cir[i], l, tmp, t)) continue;
		if (dcmp (tmp - k) < 0) k = tmp, inter = t, sign = i;
	}
	//ray 变向
	if (sign == -1) return sign;
	ray l1;
	l1.s = cir[sign].o;
	l1.dir.x = (inter.x - l1.s.x) / cir[sign].r;
	l1.dir.y = (inter.y - l1.s.y) / cir[sign].r;
	l1.dir.z = (inter.z - l1.s.z) / cir[sign].r;
	newRay (l, l1, inter);
	return sign;
}
int main ()
{
//  freopen ("in", "r", stdin);
	int i;
	scanf ("%d", &n);
	for (i = 1; i <= n; i++) //输入空间的球位置
		scanf ("%lf%lf%lf%lf", &cir[i].o.x, &cir[i].o.y, &cir[i].o.z, &cir[i].r);
	scanf ("%lf%lf%lf%lf%lf%lf", &l.s.x, &l.s.y, &l.s.z, &l.dir.x, &l.dir.y, &l.dir.z);
	for (i = 0; i <= 10; i++){ //最多输出十次相交的球的编号
		int sign = update ();
		if (sign == -1) break;
		if (i == 0) printf ("%d", sign);
		else if (i < 10) printf (" %d", sign);
		else printf (" etc.");
	}
	puts ("");
}
