#include<bits/stdc++.h>
#define MAXN 1000//����������
#define offset 10000//����������
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
#define _sign(x) ((x)>eps?1:((x)<-eps?2:0))
struct point{double x,y;};//��
struct line{point a,b;};//��
//�����ڰ�ƽ�潻
double xmult(point p1,point p2,point p0){
	return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
int same_side(point p1,point p2,point l1,point l2){
	return xmult(l1,p1,l2)*xmult(l1,p2,l2)>eps;
}
point intersection(point u1,point u2,point v1,point v2){
	point ret=u1;
	double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
			/((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
	ret.x+=(u2.x-u1.x)*t;
	ret.y+=(u2.y-u1.y)*t;
	return ret;
}
//���������l1,l2ȷ����ֱ���и���side���и�,��֤l1,l2,side������
void polygon_cut(int& n,point* p,point l1,point l2,point side){
	point pp[MAXN];
	int m=0,i;
	for (i=0;i<n;i++){
		if (same_side(p[i],side,l1,l2))
			pp[m++]=p[i];
		if (!same_side(p[i],p[(i+1)%n],l1,l2)&&!(zero(xmult(p[i],l1,l2))&&zero(xmult(p[(i+1)%n],l1,l2))))
			pp[m++]=intersection(p[i],p[(i+1)%n],l1,l2);
	}
	for (n=i=0;i<m;i++)
		if (!i||!zero(pp[i].x-pp[i-1].x)||!zero(pp[i].y-pp[i-1].y))
			p[n++]=pp[i];
	if (zero(p[n-1].x-p[0].x)&&zero(p[n-1].y-p[0].y))
		n--;
	if (n<3)
		n=0;
}
