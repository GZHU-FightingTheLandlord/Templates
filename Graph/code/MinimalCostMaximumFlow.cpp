#include <string.h>
#include <algorithm>
#include <queue>
using namespace std;

const int MAX = 1000;
const int infi = (1 << 31) - 1;

int cnt;

struct edge{
    int v, w, c, next;
    edge(){}
    edge(int vv, int ww, int cc, int nn):v(vv), w(ww), c(cc), next(nn){}
};

edge e[MAX * MAX]; // 边集
int head[MAX];     // u为起点的第一条边在边集中的位置
int pre[MAX];      // 父边
int cost[MAX];     // 费用
int inq[MAX];

void init() // 初始化
{
    cnt = 0;
    memset(head, -1, sizeof head);
}

void addedge(int u, int v, int w, int c) //加边
{
    // 正向边
    e[cnt] = edge(v, w, c, head[u]); head[u] = cnt++;
    // 反向边
    e[cnt] = edge(u, 0, -c, head[v]); head[v] = cnt++;
    // 对任意边i，i^1即为其逆向边
}

bool spfa(int from, int to)// spfa求最小费用增广路
{
    queue<int> q;
    for (int i = 0; i < MAX; i++)
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

int augment(int from, int to) // 增广
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
    while (spfa(from, to)) // 持续增广直至无法找到增广路
        ans += augment(from, to);
    return ans;
}


