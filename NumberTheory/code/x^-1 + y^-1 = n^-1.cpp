/* 
 * Source : https://www.nowcoder.com/acm/contest/90/F
 * Author : quailty
 */
#include<bits/stdc++.h>
using namespace std;
int main()
{
    int T;
    scanf("%d",&T);
    while(T--)
    {
        int n,res=1;
        scanf("%d",&n);
        for(int i=2;i*i<=n;i++)if(n%i==0)
        {
            int cnt=0;
            while(n%i==0)n/=i,cnt++;
            res*=2*cnt+1;
        }
        if(n>1)res*=3;
        printf("%d\n",(res+1)/2);
    }
    return 0;
}
