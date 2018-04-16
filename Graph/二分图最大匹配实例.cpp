/*
    Author: SemonChan
    Time: 2018-04-17 00:23
    Problem: Hihocoder-1122
    Solution Source: https://cn.vjudge.net/solution/13191150
*/
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
using namespace std;

int n, m;
vector<int> v[1002];
int link[1002];
int vis[1002];

int find(int k)
{
	for (int i = 0; i < v[k].size(); i++)
	{
		int to = v[k][i];
		if (!vis[to])
		{
			vis[to] = 1;
			if (link[to] == -1 || find(link[to]))
			{
				link[to] = k;
				link[k] = to;
				return 1;
			}
		}
	}
	return 0;
}

int solve()
{
	int ans = 0;
	memset(link, -1, sizeof link);
	for (int i = 1; i <= n;i++)
		if (link[i] == -1)
		{
			memset(vis, 0, sizeof vis);
			if (find(i)) ans++;
		}
	return ans;
}

int main()
{
	scanf("%d%d", &n, &m);
	for (int i = 1; i <= m; i++)
	{
		int a, b;
		scanf("%d%d", &a, &b);
		v[a].push_back(b);
		v[b].push_back(a);
	}
	printf("%d\n", solve());
	return 0;
}
