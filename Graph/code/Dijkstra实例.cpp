/*
    Author: SemonChan
    Time: 2018-04-17 00:06
    Problem: HDU-1874
    Solution Source: https://cn.vjudge.net/solution/13521262
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
bool vis[200];

void init()
{
    for (int i = 0; i < 200; i++)
        e[i].clear(), dis[i] = 1e9, vis[i] = 0;
}

int dijkstra(int s, int t)
{
    priority_queue<node> q;

    dis[s] = 0;
    q.push(node(s, 0));
    while (!q.empty())
    {
        int u = q.top().v;
        q.pop();

        if (u == t)
            return 1;
        if (vis[u])
            continue;
        vis[u] = 1;
        for (int i = 0; i < e[u].size(); i++)
        {
            int v = e[u][i].v, w = e[u][i].w;
            if (!vis[v] && dis[v] > dis[u] + w)
            {
                dis[v] = dis[u] + w;
                q.push(node(v, dis[v]));
            }
        }
    }
    return -1;
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
        printf("%d\n", dijkstra(s, t) == -1 ? -1 : dis[t]);
    }
    return 0;
}
