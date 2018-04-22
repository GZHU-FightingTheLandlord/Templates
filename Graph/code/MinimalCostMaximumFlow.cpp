#include <string.h>
#include <algorithm>
#include <queue>
using namespace std;

const int MAX = 1000;
const int infi = (1 << 31) - 1;

int n, m, cnt;

struct edge{
    int v, w, c, next;
    edge(){}
    edge(int vv, int ww, int cc, int nn):v(vv), w(ww), c(cc), next(nn){}
};

edge e[MAX * MAX];
int head[MAX];
int pre[MAX];
int cost[MAX];
int inq[MAX];

void init()
{
    cnt = 0;
    memset(head, -1, sizeof head);
}

void addedge(int u, int v, int w, int c)
{
    e[cnt] = edge(v, w, c, head[u]); head[u] = cnt++;
    e[cnt] = edge(u, 0, -c, head[v]); head[v] = cnt++;
}

bool spfa(int from, int to)
{
    queue<int> q;
    for (int i = 0; i < n + 5; i++)
        cost[i] = infi, pre[i] = -1, inq[i] = 0;
    q.push(from);
    inq[from] = 1; cost[from] = 0;
    while (!q.empty())
    {
        int u = q.front();
        q.pop(); inq[u] = 0;

        for (int i = head[u]; i != -1; i = e[i].next)
        {
            int v = e[i].v;
            if (e[i].w > 0 && cost[v] > cost[u] + e[i].c)
            {
                cost[v] = cost[u] + e[i].c;
                pre[v] = i;
                if (!inq[v])
                    q.push(v), inq[v] = 1;
            }
        }
    }
    return cost[to] != infi;
}

int augment(int from, int to)
{
    int now = pre[to], flow = infi;
    while (now != -1)
    {
        flow = min(flow, e[now].w);
        now = pre[e[now ^ 1].v];
    }
    now = pre[to];
    while (now != -1)
    {
        e[now].w -= flow;
        e[now ^ 1].w += flow;
        now = pre[e[now ^ 1].v];
    }
    return flow * cost[to];
}

int MinCostFlow(int from, int to)
{
    int ans = 0;
    while (spfa(from, to))
        ans += augment(from, to);
    return ans;
}


