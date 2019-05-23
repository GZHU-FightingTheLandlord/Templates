#include math.h
struct point{double x,y;};
//计算cross product (P1-P0)x(P2-P0)
double xmult(point p1,point p2,point p0){
	return (p1.x-p0.x)(p2.y-p0.y)-(p2.x-p0.x)(p1.y-p0.y);
}
//计算三角形面积,输入三顶点
double area_triangle(point p1,point p2,point p3){
	return fabs(xmult(p1,p2,p3))/2;
}
//计算三角形面积,输入三边长
double area_triangle(double a,double b,double c){
	double s=(a+b+c)/2;
	return sqrt(s(s-a)(s-b)(s-c));
}
//计算多边形面积(凸包),顶点按顺时针或逆时针给出
double area_polygon(int n,point p){
	double s1=0,s2=0;
	int i;
	for (i=0;in;i++)
		s1+=p[(i+1)%n].yp[i].x,s2+=p[(i+1)%n].yp[(i+2)%n].x;
	return fabs(s1-s2)/2;
}
