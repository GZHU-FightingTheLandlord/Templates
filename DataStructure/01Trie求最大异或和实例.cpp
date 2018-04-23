/*
  Author: SemonChan
  Problem：HDU-5536
  Solution Source: https://cn.vjudge.net/solution/13181232
*/
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <queue>
using namespace std;

struct node {
	int num;
	int next[3];
	node() { num = 0; memset(next, 0, sizeof next); }
};

int n, cnt;
node tr[100003];
const int len = 30;

void init()
{
	cnt = 0;
	tr[0] = node();
}

void put(int temp)
{
	int now = 0;
	for (int i = 0; i <= len; i++)
	{
		int k = bool((1 << (len - i)) & temp);
		if (tr[now].next[k] == 0)
		{
			tr[now].next[k] = ++cnt;
			tr[cnt] = node();
		}
		now = tr[now].next[k];
		tr[now].num++;
	}
}

int search(int temp)
{
	int now = 0, ans = 0;
	for (int i = 0; i <= len; i++)
	{
		int k = bool((1 << (len - i)) & temp);
		if (k)
		{
			if (tr[now].next[0] && tr[tr[now].next[0]].num > 0)
				now = tr[now].next[0], ans <<= 1;
			else
				now = tr[now].next[1], ans = ans << 1 | 1;
		}
		else
		{
			if (tr[now].next[1] && tr[tr[now].next[1]].num > 0)
				now = tr[now].next[1], ans = ans << 1 | 1;
			else
				now = tr[now].next[0], ans <<= 1;
		}
	}
	return temp ^ ans;
}

void del(int temp)
{
	int now = 0;
	for (int i = 0; i <= len; i++)
	{
		now = tr[now].next[bool((1 << (len - i)) & temp)];
		tr[now].num--;
	}
}

int s[1010];

int main()
{
	int t, n;
	scanf("%d", &t);
	while (t--)
	{
		init();
		scanf("%d", &n);
		for (int i = 1; i <= n; i++)
			scanf("%d", &s[i]), put(s[i]);
		int ma = -1111;
		for (int i = 1; i < n; i++)
		{
			del(s[i]);
			for (int j = i + 1; j <= n; j++)
			{
				del(s[j]);
				ma = max(ma, search(s[i] + s[j]));
				put(s[j]);
			}
			put(s[i]);
		}
		printf("%d\n", ma);
	}
	return 0;
}
