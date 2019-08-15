const int N=1e5+7;
struct point{
    double x,y;
}S[N];
inline double Cross(point a,point b,point c){
    return (a.x-c.x)*(b.y-c.y)-(a.y-c.y)*(b.x-c.x);
}
inline double Dis(point a,point b){
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
//先求凸包再旋转卡壳
double FarthestPointPair(int top){
    if(top==1)return Dis(S[0],S[1]); 
    S[++top]=S[0];int j=2;double ans=0;
    for(int i=0;i<top;++i){
        while(Cross(S[i],S[i+1],S[j])<Cross(S[i],S[i+1],S[j+1]))
            j=(j+1)%top;
        ans=max(ans,max(Dis(S[i],S[j]),Dis(S[i+1],S[j])));
    }
    return ans;
}