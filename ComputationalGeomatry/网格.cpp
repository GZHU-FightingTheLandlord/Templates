#define abs(x) ((x)>0?(x):-(x))
struct point{int x,y;};
int gcd(int a,int b){return b?gcd(b,a%b):a;}
//������ϵ���������
int grid_onedge(int n,point* p){
	int i,ret=0;
	for (i=0;i<n;i++)
		ret+=gcd(abs(p[i].x-p[(i+1)%n].x),abs(p[i].y-p[(i+1)%n].y));
	return ret;
}
//������ڵ���������
int grid_inside(int n,point* p){
	int i,ret=0;
	for (i=0;i<n;i++)
		ret+=p[(i+1)%n].y*(p[i].x-p[(i+2)%n].x);
	return (abs(ret)-grid_onedge(n,p))/2+1;
}
