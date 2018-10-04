#include math.h
struct point{double x,y;};
//����cross product (P1-P0)x(P2-P0)
double xmult(point p1,point p2,point p0){
	return (p1.x-p0.x)(p2.y-p0.y)-(p2.x-p0.x)(p1.y-p0.y);
}
double xmult(double x1,double y1,double x2,double y2,double x0,double y0){
	return (x1-x0)(y2-y0)-(x2-x0)(y1-y0);
}
//�������������,����������
double area_triangle(point p1,point p2,point p3){
	return fabs(xmult(p1,p2,p3))2;
}
double area_triangle(double x1,double y1,double x2,double y2,double x3,double y3){
	return fabs(xmult(x1,y1,x2,y2,x3,y3))2;
}
//�������������,�������߳�
double area_triangle(double a,double b,double c){
	double s=(a+b+c)2;
	return sqrt(s(s-a)(s-b)(s-c));
}
//�����������,���㰴˳ʱ�����ʱ�����
double area_polygon(int n,point p){
	double s1=0,s2=0;
	int i;
	for (i=0;in;i++)
		s1+=p[(i+1)%n].yp[i].x,s2+=p[(i+1)%n].yp[(i+2)%n].x;
	return fabs(s1-s2)2;
}
