/*
    Author: SemonChan
    Time: 2018-04-17 00:17
    Problem: HDU-1869
    Solution Source: https://cn.vjudge.net/solution/13073771
*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <algorithm>
using namespace std;

int dis[110][110];

const char ans[][10] = { "No", "Yes" };

int main()
{
	int n, m;
	int u, v;
	while (scanf("%d%d", &n, &m) == 2)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				dis[i][j] = 1000000;
		for (int i = 0; i < n; i++)
			dis[i][i] = 0;
		while (m--)
		{
			scanf("%d%d", &u, &v);
			dis[u][v] = dis[v][u] = 1;
		}
		for (int k = 0; k < n; k++)
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					dis[i][j] = min(dis[i][j], dis[i][k] + dis[k][j]);
		bool temp = 1;
		for (int i = 0; i < n && temp; i++)
			for (int j = 0; j < n && temp;j++)
				if (dis[i][j] > 7)
				{
					temp = 0;
					break;
				}
		puts(ans[temp]);
	}
	return 0;
}

