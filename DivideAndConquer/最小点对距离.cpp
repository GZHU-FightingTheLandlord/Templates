#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

struct node{
    double x, y;
    node(double x_ = 0, double y_ = 0):x(x_), y(y_){}
    bool operator < (const node b)const
    {
        return y < b.y;
    }
};

node v[100005];
node nd[100005];

/* 求欧氏距离    */
inline double getd(const node a, const node b)
{
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double solve(int a, int b)
{
    /* 若区间只有一个点   */
    if (a == b)
        return 1e30;

    int mid = (a + b) / 2, cnt = 0;
    double d = min(solve(a, mid), solve(mid + 1, b)); // 递归求两个子区间的最小点对距离

    /* 遍历区间，搜索与v[mid]横坐标差值小于d的点 */
    for (int i = mid; i >= a; i--)
        if (fabs(v[i].x - v[mid].x) < d)
            nd[cnt++] = v[i];
        else
            break;
    for (int i = mid + 1; i <= b; i++)
        if (fabs(v[i].x - v[mid].x) < d)
            nd[cnt++] = v[i];
        else
            break;

    /* 按纵坐标从小到大排序   */
    sort(nd, nd + cnt);

    /* 枚举求区间内最小点对距离 */
    for (int i = 0; i < cnt - 1; i++)
        for (int j = i + 1; j < cnt; j++)
        {
            if (fabs(nd[i].y - nd[j].y) > d)
                break;
            d = min(d, getd(nd[i], nd[j]));
        }
    return d;
}

bool cmp(const node a, const node b)
{
    return a.x < b.x;
}

int main()
{
    int n;
    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
        scanf("%lf%lf", &v[i].x, &v[i].y);
    
    /* 按横坐标从小到大排序   */
    sort(v + 1, v + 1 + n, cmp);
    printf("%.2f\n", solve(1, n));
    return 0;
}
