const double eps = 1e-8;            //浮点数精度控制
struct point                        //点或者向量结构
{
    double x,y;
    point(double _x=0.0,double _y=0.0)
        : x(_x),y(_y) {}
    point operator - (const point & v)
    {
        return point(x-v.x,y-v.y);
    }
    double sqrx()                    //向量的模
    {
        return sqrt(x*x+y*y);
    }
};
double xmult(point & p1,point & p2,point & p0)        //叉乘
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x);
}
double distancex(point & p1,point & p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
point intersection(point u1,point u2,point v1,point v2)        //两直线交点
{
    point ret=u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
            /((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x+=(u2.x-u1.x)*t;
    ret.y+=(u2.y-u1.y)*t;
    return ret;
}
void intersection_line_circle(point c,double r,point l1,point l2,point& p1,point& p2){
    point p=c;
    double t;
    p.x+=l1.y-l2.y;
    p.y+=l2.x-l1.x;
    p=intersection(p,c,l1,l2);
    t=sqrt(r*r-distancex(p,c)*distancex(p,c))/distancex(l1,l2);
    p1.x=p.x+(l2.x-l1.x)*t;
    p1.y=p.y+(l2.y-l1.y)*t;
    p2.x=p.x-(l2.x-l1.x)*t;
    p2.y=p.y-(l2.y-l1.y)*t;
}
point ptoseg(point p,point l1,point l2)            //点到线段的最近距离
{
    point t=p;
    t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
    if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
    return distancex(p,l1)<distancex(p,l2)?l1:l2;
    return intersection(p,t,l1,l2);
}
double distp(point & a,point & b)
{
    return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y);
}
double Direct_Triangle_Circle_Area(point a,point b,point o,double r)
{
    double sign=1.0;
    a=a-o;
    b=b-o;
    o=point(0.0,0.0);
    if(fabs(xmult(a,b,o))<eps) return 0.0;
    if(distp(a,o)>distp(b,o))
    {
        swap(a,b);
        sign=-1.0;
    }
    if(distp(a,o)<r*r+eps)
    {
        if(distp(b,o)<r*r+eps) return xmult(a,b,o)/2.0*sign;
        point p1,p2;
        intersection_line_circle(o,r,a,b,p1,p2);
        if(distancex(p1,b)>distancex(p2,b)) swap(p1,p2);
        double ret1=fabs(xmult(a,p1,o));
        double ret2=acos( p1*b/p1.sqrx()/b.sqrx() )*r*r;
        double ret=(ret1+ret2)/2.0;
        if(xmult(a,b,o)<eps && sign>0.0 || xmult(a,b,o)>eps && sign<0.0) ret=-ret;
        return ret;
    }
    point ins=ptoseg(o,a,b);
    if(distp(o,ins)>r*r-eps)
    {
        double ret=acos( a*b/a.sqrx()/b.sqrx() )*r*r/2.0;
        if(xmult(a,b,o)<eps && sign>0.0 || xmult(a,b,o)>eps && sign<0.0) ret=-ret;
        return ret;
    }
    point p1,p2;
    intersection_line_circle(o,r,a,b,p1,p2);
    double cm=r/(distancex(o,a)-r);
    point m=point( (o.x+cm*a.x)/(1+cm) , (o.y+cm*a.y)/(1+cm) );
    double cn=r/(distancex(o,b)-r);
    point n=point( (o.x+cn*b.x)/(1+cn) , (o.y+cn*b.y)/(1+cn) );
    double ret1 = acos( m*n/m.sqrx()/n.sqrx() )*r*r;
    double ret2 = acos( p1*p2/p1.sqrx()/p2.sqrx() )*r*r-fabs(xmult(p1,p2,o));
    double ret=(ret1-ret2)/2.0;
    if(xmult(a,b,o)<eps && sign>0.0 || xmult(a,b,o)>eps && sign<0.0) ret=-ret;
    return ret;
}
