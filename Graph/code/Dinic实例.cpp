/*
    Author: SemonChan
    Time: 2018-04-17 00:26
    Problem: POJ-3281
    Solution Source: https://cn.vjudge.net/solution/12677665
*/
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <queue>
using namespace std;

const int MAX = 420;

struct ee {
	int t, c;
	int next;
}edge[MAX * 100];

int n, cnt;
int head[MAX];
int level[MAX];
int cur[MAX];

void addedge(int f, int t)
{
	edge[cnt].t = t;
	edge[cnt].c = 1;
	edge[cnt].next = head[f];
	head[f] = cnt++;

	edge[cnt].t = f;
	edge[cnt].c = 0;
	edge[cnt].next = head[t];
	head[t] = cnt++;
}

bool bfs()
{
	queue<int> q;
	while (!q.empty())
		q.pop();
	memset(level, 0, sizeof level);

	q.push(0);
	level[0] = 1;
	while (!q.empty())
	{
		int t = q.front();
		q.pop();
		for (int i = head[t]; i != -1; i = edge[i].next)
		{
			int tt = edge[i].t;
			if (!level[tt] && edge[i].c)
			{
				q.push(tt);
				level[tt] = level[t] + 1;
			}
		}
	}
	return level[n] != 0;
}

int dfs(int t)
{
	if (t == n)
		return 1;
	for (int& i = cur[t]; i != -1; i = edge[i].next)
	{
		struct ee& e = edge[i];
		if (level[e.t] == level[t] + 1 && e.c && dfs(e.t))
		{
			e.c = 0;
			edge[i ^ 1].c = 1;
			return 1;
		}
	}
	return 0;
}

int dinic()
{
	int ans = 0;
	while (bfs())
	{
		for (int i = 0; i <= n; i++)
			cur[i] = head[i];
		while (dfs(0))
			ans++;
	}
	return ans;
}

int main()
{
	int N, f, t;
	while (scanf("%d%d%d", &N, &f, &t) == 3)
	{
		n = f + t + 2 * N + 1;
		cnt = 0;
		for (int i = 0; i <= n; i++)
			head[i] = -1;
		for (int i = 1; i <= f; i++)
			addedge(0, i);
		for (int i = 1; i <= t; i++)
			addedge(f + 2 * N + i, n);
		for (int i = 1; i <= N; i++)
			addedge(f + 2 * i - 1, f + 2 * i);
		
		for (int i = 1; i <= N; i++)
		{
			int ff, tt;
			scanf("%d%d", &ff, &tt);
			while (ff--)
			{
				int k;
				scanf("%d", &k);
				addedge(k, f + 2 * i - 1);
			}
			while (tt--)
			{
				int k;
				scanf("%d", &k);
				addedge(f + 2 * i, f + 2 * N + k);
			}
		}
		printf("%d\n", dinic());
	}
	return 0;
}
