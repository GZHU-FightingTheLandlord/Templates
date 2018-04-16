/*
    Author: SemonChan
    Time: 2018-04-17 00:15
    Problem: HDU-1874
    Solution Source: https://cn.vjudge.net/solution/13521322
*/

#include <bits/stdc++.h>
using namespace std;

struct node{
    int v, w;
    node(int v_ = 0, int w_ = 0) : v(v_), w(w_){}
    bool operator< (const node b)const
    {
        return w > b.w;
    }
};

int n, m;
vector<node> e[200];
int dis[200];
bool inq[200];

void init()
{
    for (int i = 0; i < 200; i++)
        e[i].clear(), dis[i] = 1e9, inq[i] = 0;
}

int spfa(int s, int t)
{
    queue<int> q;

    q.push(s);
    dis[s] = 0; inq[s] = 1;
    while (!q.empty())
    {
        int u = q.front();
        q.pop(); inq[u] = 0;

        for (int i = 0; i < e[u].size(); i++)
        {
            int v = e[u][i].v, w = e[u][i].w;
            if (dis[v] > dis[u] + w)
            {
                dis[v] = dis[u] + w;
                if (!inq[v])
                    q.push(v);
            }
        }
    }
    return dis[t] == 1e9 ? -1 : 1;
}

int main()
{
    while (scanf("%d%d", &n, &m) == 2)
    {
        init();
        while (m--)
        {
            int u, v, w;
            scanf("%d%d%d", &u, &v, &w);
            e[u].push_back(node(v, w));
            e[v].push_back(node(u, w));
        }
        int s, t;
        scanf("%d%d", &s, &t);
        printf("%d\n", spfa(s, t) == -1 ? -1 : dis[t]);
    }
    return 0;
}
