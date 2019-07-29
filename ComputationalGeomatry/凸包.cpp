const int N=10005;
//a[N]:待处理点,n:待处理点数,p[N]:凸包点,top:凸包点数
int top,n;
struct node{
    int x,y;
}a[N],p[N];
double cross(node p0,node p1,node p2){
    return (p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x);
}
double dis(node a,node b){
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
//z=0代表三点在一条直线上,a[0]代表旋转中心点
bool cmp(node p1,node p2){
    double z=cross(a[0],p1,p2);
    if(z>0||(z==0&&dis(a[0],p1)<dis(a[0],p2)))
        return true;
    return false;
}
void Graham(){
    int k=0;
    for(int i=0; i<n; ++i)
        if(a[i].y<a[k].y||(a[i].y==a[k].y&&a[i].x<a[k].x))
            k=i;
    swap(a[0],a[k]);
    sort(a+1,a+n,cmp);
    top=1,p[0]=a[0],p[1]=a[1];
    for(int i=2; i<n; ++i){
        while(cross(p[top-1],p[top],a[i])<0&&top)
            top--;
        p[++top]=a[i];
    }
}
