#include<bits/stdc++.h>
#define MAXN 1000//����������
#define offset 10000//����������
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
#define _sign(x) ((x)>eps?1:((x)<-eps?2:0))
struct point{double x,y;};//��
struct line{point a,b;};//��
//���
double xmult(point p1,point p2,point p0){
	return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
//�ж�͹�����,���㰴˳ʱ�����ʱ�����,�������ڱ߹���
int is_convex(int n,point* p){
	int i,s[3]={1,1,1};
	for (i=0;i<n&&s[1]|s[2];i++)
		s[_sign(xmult(p[(i+1)%n],p[(i+2)%n],p[i]))]=0;
	return s[1]|s[2];
}
//�ж�͹�����,���㰴˳ʱ�����ʱ�����,���������ڱ߹���
int is_convex_v2(int n,point* p){
	int i,s[3]={1,1,1};
	for (i=0;i<n&&s[0]&&s[1]|s[2];i++)
		s[_sign(xmult(p[(i+1)%n],p[(i+2)%n],p[i]))]=0;
	return s[0]&&s[1]|s[2];
}
//�е���͹������ڻ����α���,���㰴˳ʱ�����ʱ�����
int inside_convex(point q,int n,point* p){
	int i,s[3]={1,1,1};
	for (i=0;i<n&&s[1]|s[2];i++)
		s[_sign(xmult(p[(i+1)%n],q,p[i]))]=0;
	return s[1]|s[2];
}
//�е���͹�������,���㰴˳ʱ�����ʱ�����,�ڶ���α��Ϸ���0
int inside_convex_v2(point q,int n,point* p){
	int i,s[3]={1,1,1};
	for (i=0;i<n&&s[0]&&s[1]|s[2];i++)
		s[_sign(xmult(p[(i+1)%n],q,p[i]))]=0;
	return s[0]&&s[1]|s[2];
}
//�е�������������,���㰴˳ʱ�����ʱ�����
//on_edge��ʾ���ڶ���α���ʱ�ķ���ֵ
int inside_polygon(point q,int n,point* p,int on_edge=1){
	point q2;
	int i=0,count;
	while (i<n)
		for (count=i=0,q2.x=rand()+offset,q2.y=rand()+offset;i<n;i++)
			if (zero(xmult(q,p[i],p[(i+1)%n]))&&(p[i].x-q.x)*(p[(i+1)%n].x-q.x)<eps&&(p[i].y-q.y)*(p[(i+1)%n].y-q.y)<eps)
				return on_edge;
			else if (zero(xmult(q,q2,p[i])))
				break;
			else if (xmult(q,p[i],q2)*xmult(q,p[(i+1)%n],q2)<-eps&&xmult(p[i],q,p[(i+1)%n])*xmult(p[i],q2,p[(i+1)%n])<-eps)
				count++;
	return count&1;
}
inline int opposite_side(point p1,point p2,point l1,point l2){
	return xmult(l1,p1,l2)*xmult(l1,p2,l2)<-eps;
}

inline int dot_online_in(point p,point l1,point l2){
	return zero(xmult(p,l1,l2))&&(l1.x-p.x)*(l2.x-p.x)<eps&&(l1.y-p.y)*(l2.y-p.y)<eps;
}
//���߶�������������,���㰴˳ʱ�����ʱ�����,��߽��ཻ����1
int inside_polygon(point l1,point l2,int n,point* p){
	point t[MAXN],tt;
	int i,j,k=0;
	if (!inside_polygon(l1,n,p)||!inside_polygon(l2,n,p))
		return 0;
	for (i=0;i<n;i++)
		if (opposite_side(l1,l2,p[i],p[(i+1)%n])&&opposite_side(p[i],p[(i+1)%n],l1,l2))
			return 0;
		else if (dot_online_in(l1,p[i],p[(i+1)%n]))
			t[k++]=l1;
		else if (dot_online_in(l2,p[i],p[(i+1)%n]))
			t[k++]=l2;
		else if (dot_online_in(p[i],l1,l2))
			t[k++]=p[i];
	for (i=0;i<k;i++)
		for (j=i+1;j<k;j++){
			tt.x=(t[i].x+t[j].x)/2;
			tt.y=(t[i].y+t[j].y)/2;
			if (!inside_polygon(tt,n,p))
				return 0;
		}
	return 1;
}
point intersection(line u,line v){
	point ret=u.a;
	double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
			/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
	ret.x+=(u.b.x-u.a.x)*t;
	ret.y+=(u.b.y-u.a.y)*t;
	return ret;
}
point barycenter(point a,point b,point c){
	line u,v;
	u.a.x=(a.x+b.x)/2;
	u.a.y=(a.y+b.y)/2;
	u.b=c;
	v.a.x=(a.x+c.x)/2;
	v.a.y=(a.y+c.y)/2;
	v.b=b;
	return intersection(u,v);
}
//���������
point barycenter(int n,point* p){
	point ret,t;
	double t1=0,t2;
	int i;
	ret.x=ret.y=0;
	for (i=1;i<n-1;i++)
		if (fabs(t2=xmult(p[0],p[i],p[i+1]))>eps){
			t=barycenter(p[0],p[i],p[i+1]);
			ret.x+=t.x*t2;
			ret.y+=t.y*t2;
			t1+=t2;
		}
	if (fabs(t1)>eps)
		ret.x/=t1,ret.y/=t1;
	return ret;
}
