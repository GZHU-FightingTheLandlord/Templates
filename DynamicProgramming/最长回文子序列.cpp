/*
  Source : https://www.nowcoder.com/acm/contest/90/D
  Author : Semon999
*/

#include <bits/stdc++.h>
using namespace std;

#define pb 		push_back
#define mp 		make_pair
#define fi 		first
#define se 		second
#define all(a)  (a).begin(),(a).end()
#define szz(a) 	(int)a.size()

typedef long long ll;
typedef long double ld;
typedef double db;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef pair<ll, int> pli;
typedef pair<db, db> pdd;

int len;
char v[1300];
int dp[1300][1300];

int dfs(int l, int r)
{
	if (l == r) return 1;
	if (l > r) return 0;
	if (dp[l][r] != -1) return dp[l][r];
	int res = max(dfs(l, r - 1), dfs(l + 1, r));
	if (v[l] == v[r]) res = max(res, dfs(l + 1, r - 1) + 2);
	return dp[l][r] = res;
}

int main()
{
	while (scanf("%s", v) != EOF) {
		len = strlen(v);
		for (int i = 0; i < len; i++) v[i] = tolower(v[i]);
		memset(dp, -1, sizeof dp);
		printf("%d\n", len - dfs(0, len - 1));
	}
}
