/*
    Author: SemonChan
    Time: 2018-04-17 00:19
    Problem: HDU-1875
    Solution Source: https://cn.vjudge.net/solution/12270502
*/
#include<bits/stdc++.h>
using namespace std;

struct edge{
    int x, y;
    double l;
}e[5000];
int n, edm;
int pre[110];
double ans;
int mapp[110][4];

bool cmp(const edge a, const edge b)
{
    return a.l < b.l;
}

void ffind(int &x)
{
    while (x != pre[x])
        x = pre[x];
}

void kruskal()
{
    int i;
    sort(e + 1, e + 1 + edm, cmp);
    for (i = 1;i <= edm;i++)
    {
        ffind(e[i].x);
        ffind(e[i].y);
        if (e[i].x == e[i].y)
            continue;
        pre[e[i].y] = e[i].x;
        ans += e[i].l;
    }
}

inline double getdir(int a,int b,int c,int d)
{
    return sqrt(pow(fabs(a - c), 2) + pow(fabs(b - d), 2));
}

int main()
{
    int i,j,t;
    double dir;
    int temp;
    scanf("%d", &t);
    while (t--)
    {
        scanf("%d", &n);
        ans = 0;
        for (i = 1;i <= n;i++)
            pre[i] = i;
        for (i = 1;i <= n;i++)
            scanf("%d%d", &mapp[i][1], &mapp[i][2]);
        edm = 0;
        for (i = 1;i <= n;i++)
        {
            for (j = i + 1;j <= n;j++)
            {
                dir = getdir(mapp[i][1],mapp[i][2],mapp[j][1],mapp[j][2]);
                if (dir >= 10 && dir <= 1000)
                {
                    e[++edm].x = i;
                    e[edm].y = j;
                    e[edm].l = dir;
                }
            }
        }
        kruskal();
        temp = 0;
        for (i = 1;i <= n;i++)
            if (pre[i] == i)
                temp++;
        if (temp <= 1)
            printf("%.1f\n", ans * 100);
        else
            printf("oh!\n");
    }
    return 0;
}
