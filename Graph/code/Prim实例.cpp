/*
    Author: SemonChan
    Time: 2018-04-17 00:21
    Problem: Hihocoder-1109
    Solution Source: https://cn.vjudge.net/solution/13030273
*/
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <queue>
using namespace std;

struct node {
	int num, dis;
	node(){}
	node(int a, int b) : num(a), dis(b){}
	bool operator< (const node b)const 
	{
		return dis > b.dis;
	}
};

struct edge {
	int v, w;
	int next;
	edge(){}
};

int n, m, cnt;
long long ans;
int head[100005];
edge e[2000002];
bool vis[100005];
int dist[100005];

inline void init()
{
	cnt = ans = 0;
	memset(head, -1, (n + 2) * sizeof(int));
	memset(vis, 0, (n + 2) * sizeof(bool));
	memset(dist, 0x3f3f3f3f, (n + 2) * sizeof(int));
}

inline void addedge(int u, int v, int w)
{
	e[++cnt].v = v; e[cnt].w = w; e[cnt].next = head[u]; head[u] = cnt;
	e[++cnt].v = u; e[cnt].w = w; e[cnt].next = head[v]; head[v] = cnt;
}

void prim()
{
	priority_queue<node> q;
	while (!q.empty())
		q.pop();

	q.push(node(1, 0));
	dist[1] = 0;
	while (n)
	{
		int temp = q.top().num, dis = q.top().dis;
		q.pop();

		if (vis[temp])
			continue;
		vis[temp] = 1;
		dist[temp] = 0;
		--n;
		ans += dis;
		for (int i = head[temp]; ~i; i = e[i].next)
			if (!vis[e[i].v] && dist[e[i].v] > e[i].w)
			{
				q.push(node(e[i].v, e[i].w));
				dist[e[i].v] = e[i].w;
			}
	}
}

int main()
{
	scanf("%d%d", &n, &m);
	init();
	while (m--)
	{
		int u, v, w;
		scanf("%d%d%d", &u, &v, &w);
		addedge(u, v, w);
	}
	prim();
	printf("%lld\n", ans);
	return 0;
}


