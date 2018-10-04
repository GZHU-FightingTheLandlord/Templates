//浮点几何函数库
#include <math.h>
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
struct point{double x,y;};
struct line{point a,b;};
//计算cross product (P1-P0)x(P2-P0)
double xmult(point p1,point p2,point p0){
	return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
double xmult(double x1,double y1,double x2,double y2,double x0,double y0){
	return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
}
//计算dot product (P1-P0).(P2-P0)
double dmult(point p1,point p2,point p0){
	return (p1.x-p0.x)*(p2.x-p0.x)+(p1.y-p0.y)*(p2.y-p0.y);
}
double dmult(double x1,double y1,double x2,double y2,double x0,double y0){
	return (x1-x0)*(x2-x0)+(y1-y0)*(y2-y0);
}
//两点距离
double distance(point p1,point p2){
	return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
double distance(double x1,double y1,double x2,double y2){
	return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
//判三点共线
int dots_inline(point p1,point p2,point p3){
	return zero(xmult(p1,p2,p3));
}
int dots_inline(double x1,double y1,double x2,double y2,double x3,double y3){
	return zero(xmult(x1,y1,x2,y2,x3,y3));
}
//判点是否在线段上,包括端点
int dot_online_in(point p,line l){
	return zero(xmult(p,l.a,l.b))&&(l.a.x-p.x)*(l.b.x-p.x)<eps&&(l.a.y-p.y)*(l.b.y-p.y)<eps;
}
int dot_online_in(point p,point l1,point l2){
	return zero(xmult(p,l1,l2))&&(l1.x-p.x)*(l2.x-p.x)<eps&&(l1.y-p.y)*(l2.y-p.y)<eps;
}
int dot_online_in(double x,double y,double x1,double y1,double x2,double y2){
	return zero(xmult(x,y,x1,y1,x2,y2))&&(x1-x)*(x2-x)<eps&&(y1-y)*(y2-y)<eps;
}
//判点是否在线段上,不包括端点
int dot_online_ex(point p,line l){
	return dot_online_in(p,l)&&(!zero(p.x-l.a.x)||!zero(p.y-l.a.y))&&(!zero(p.x-l.b.x)||!zero(p.y-l.b.y));
}
int dot_online_ex(point p,point l1,point l2){
	return dot_online_in(p,l1,l2)&&(!zero(p.x-l1.x)||!zero(p.y-l1.y))&&(!zero(p.x-l2.x)||!zero(p.y-l2.y));
}
int dot_online_ex(double x,double y,double x1,double y1,double x2,double y2){
	return dot_online_in(x,y,x1,y1,x2,y2)&&(!zero(x-x1)||!zero(y-y1))&&(!zero(x-x2)||!zero(y-y2));
}
//判两点在线段同侧,点在线段上返回0
int same_side(point p1,point p2,line l){
	return xmult(l.a,p1,l.b)*xmult(l.a,p2,l.b)>eps;
}
int same_side(point p1,point p2,point l1,point l2){
	return xmult(l1,p1,l2)*xmult(l1,p2,l2)>eps;
}
//判两点在线段异侧,点在线段上返回0
int opposite_side(point p1,point p2,line l){
	return xmult(l.a,p1,l.b)*xmult(l.a,p2,l.b)<-eps;
}
int opposite_side(point p1,point p2,point l1,point l2){
	return xmult(l1,p1,l2)*xmult(l1,p2,l2)<-eps;
}
//判两直线平行
int parallel(line u,line v){
	return zero((u.a.x-u.b.x)*(v.a.y-v.b.y)-(v.a.x-v.b.x)*(u.a.y-u.b.y));
}
int parallel(point u1,point u2,point v1,point v2){
	return zero((u1.x-u2.x)*(v1.y-v2.y)-(v1.x-v2.x)*(u1.y-u2.y));
}
//判两直线垂直
int perpendicular(line u,line v){
	return zero((u.a.x-u.b.x)*(v.a.x-v.b.x)+(u.a.y-u.b.y)*(v.a.y-v.b.y));
}
int perpendicular(point u1,point u2,point v1,point v2){
	return zero((u1.x-u2.x)*(v1.x-v2.x)+(u1.y-u2.y)*(v1.y-v2.y));
}
//判两线段相交,包括端点和部分重合
int intersect_in(line u,line v){
	if (!dots_inline(u.a,u.b,v.a)||!dots_inline(u.a,u.b,v.b))
		return !same_side(u.a,u.b,v)&&!same_side(v.a,v.b,u);
	return dot_online_in(u.a,v)||dot_online_in(u.b,v)||dot_online_in(v.a,u)||dot_online_in(v.b,u);
}
int intersect_in(point u1,point u2,point v1,point v2){
	if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
		return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2);
	return dot_online_in(u1,v1,v2)||dot_online_in(u2,v1,v2)||dot_online_in(v1,u1,u2)||dot_online_in(v2,u1,u2);
}
//判两线段相交,不包括端点和部分重合
int intersect_ex(line u,line v){
	return opposite_side(u.a,u.b,v)&&opposite_side(v.a,v.b,u);
}
int intersect_ex(point u1,point u2,point v1,point v2){
	return opposite_side(u1,u2,v1,v2)&&opposite_side(v1,v2,u1,u2);
}
//计算两直线交点,注意事先判断直线是否平行!
//线段交点请另外判线段相交(同时还是要判断是否平行!)
point intersection(line u,line v){
	point ret=u.a;
	double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
			/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
	ret.x+=(u.b.x-u.a.x)*t;
	ret.y+=(u.b.y-u.a.y)*t;
	return ret;
}
point intersection(point u1,point u2,point v1,point v2){
	point ret=u1;
	double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
			/((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
	ret.x+=(u2.x-u1.x)*t;
	ret.y+=(u2.y-u1.y)*t;
	return ret;
}
//点到直线上的最近点
point ptoline(point p,line l){
	point t=p;
	t.x+=l.a.y-l.b.y,t.y+=l.b.x-l.a.x;
	return intersection(p,t,l.a,l.b);
}
point ptoline(point p,point l1,point l2){
	point t=p;
	t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
	return intersection(p,t,l1,l2);
}
//点到直线距离
double disptoline(point p,line l){
	return fabs(xmult(p,l.a,l.b))/distance(l.a,l.b);
}
double disptoline(point p,point l1,point l2){
	return fabs(xmult(p,l1,l2))/distance(l1,l2);
}
double disptoline(double x,double y,double x1,double y1,double x2,double y2){
	return fabs(xmult(x,y,x1,y1,x2,y2))/distance(x1,y1,x2,y2);
}
//点到线段上的最近点
point ptoseg(point p,line l){
	point t=p;
	t.x+=l.a.y-l.b.y,t.y+=l.b.x-l.a.x;
	if (xmult(l.a,t,p)*xmult(l.b,t,p)>eps)
		return distance(p,l.a)<distance(p,l.b)?l.a:l.b;
	return intersection(p,t,l.a,l.b);
}
point ptoseg(point p,point l1,point l2){
	point t=p;
	t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
	if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
		return distance(p,l1)<distance(p,l2)?l1:l2;
	return intersection(p,t,l1,l2);
}
//点到线段距离
double disptoseg(point p,line l){
	point t=p;
	t.x+=l.a.y-l.b.y,t.y+=l.b.x-l.a.x;
	if (xmult(l.a,t,p)*xmult(l.b,t,p)>eps)
		return distance(p,l.a)<distance(p,l.b)?distance(p,l.a):distance(p,l.b);
	return fabs(xmult(p,l.a,l.b))/distance(l.a,l.b);
}
double disptoseg(point p,point l1,point l2){
	point t=p;
	t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
	if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
		return distance(p,l1)<distance(p,l2)?distance(p,l1):distance(p,l2);
	return fabs(xmult(p,l1,l2))/distance(l1,l2);
}
//矢量V以P为顶点逆时针旋转angle并放大scale倍
point rotate(point v,point p,double angle,double scale){
	point ret=p;
	v.x-=p.x,v.y-=p.y;
	p.x=scale*cos(angle);
	p.y=scale*sin(angle);
	ret.x+=v.x*p.x-v.y*p.y;
	ret.y+=v.x*p.y+v.y*p.x;
	return ret;
}
//p点关于直线L的对称点
ponit symmetricalPointofLine(point p, line L)
{
    point p2;
    double d;
    d = L.a * L.a + L.b * L.b;
    p2.x = (L.b * L.b * p.x - L.a * L.a * p.x -
            2 * L.a * L.b * p.y - 2 * L.a * L.c) / d;
    p2.y = (L.a * L.a * p.y - L.b * L.b * p.y -
            2 * L.a * L.b * p.x - 2 * L.b * L.c) / d;
    return p2;
}
//求两点的平分线
line bisector(point& a, point& b) {
	line ab, ans;  ab.set(a, b);
	double midx = (a.x + b.x)/2.0,	midy = (a.y + b.y)/2.0;
	ans.a = -ab.b, ans.b = -ab.a, ans.c = -ab.b * midx + ab.a * midy;
	return ans;
}
// 已知入射线、镜面，求反射线。
// a1,b1,c1为镜面直线方程(a1 x + b1 y + c1 = 0 ,下同)系数;
a2,b2,c2为入射光直线方程系数;
a,b,c为反射光直线方程系数.
// 光是有方向的，使用时注意：入射光向量:<-b2,a2>；反射光向量:<b,-a>.
// 不要忘记结果中可能会有"negative zeros"
void reflect(double a1,double b1,double c1,
double a2,double b2,double c2,
double &a,double &b,double &c)
{
	double n,m;
	double tpb,tpa;
	tpb=b1*b2+a1*a2;
	tpa=a2*b1-a1*b2;
	m=(tpb*b1+tpa*a1)/(b1*b1+a1*a1);
	n=(tpa*b1-tpb*a1)/(b1*b1+a1*a1);
	if(fabs(a1*b2-a2*b1)<1e-20)
	{
		a=a2;b=b2;c=c2;
		return;
	}
	double xx,yy; //(xx,yy)是入射线与镜面的交点。
	xx=(b1*c2-b2*c1)/(a1*b2-a2*b1);
	yy=(a2*c1-a1*c2)/(a1*b2-a2*b1);
	a=n;
	b=-m;
	c=m*yy-xx*n;
}
