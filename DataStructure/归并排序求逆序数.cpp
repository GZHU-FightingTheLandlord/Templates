/*
	Author: Seee
	Problem Source: POJ-2299
	Code Source: https://cn.vjudge.net/solution/13942229
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <bitset>
#include <math.h>
#include <cmath>
#include <time.h>
#include <string.h>
#include <vector>
#include <set>
#include <deque>
#include <stack>
#include <time.h>
#include <map>
#include <queue>
#include <functional>
#include <cctype>
#include <iomanip>

#pragma warning(disable:4996)
#pragma comment(linker, "/STACK:336777216")
using namespace std;

#define pb(a)        push_back(a)
#define mp(a, b)     make_pair(a, b)
#define all(a)       a.begin(),a.end()
#define szz(a)       (int)a.size()
#define endl         '\n'

typedef long long ll;
typedef double db;
typedef long double ld;
typedef unsigned long long ull;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
typedef pair<ll, int> pli;
typedef pair<db, db> pdd;

const ll MOD = 1e9 + 7;
const int inf = 0x3f3f3f3f;
const int INF = 0x7fffffff;
const ll llINF = 0x3f3f3f3f3f3f3f3f;
const ll LLINF = 0x7fffffffffffffff;
const db PI = acos(-1.0);
const db eps = 1e-8;

inline ll g() { char c; ll p = 1; while (!isdigit(c = getchar())) if (c == '-') p = -p; ll x = c - '0'; while (isdigit(c = getchar())) x = x * 10 + c - '0'; return p * x; }

int v[600005];
int tmp[600005];
ll ans;

void merge(int l, int mid, int r)
{
	int num = mid - l + 1;
	int i = l, j = mid + 1, k = l;
	while (i <= mid && j <= r) {
		if (v[i] <= v[j]) {
			tmp[k++] = v[i++];
			num--;
		}
		else {
			tmp[k++] = v[j++];
			ans += num;
		}
	}
	while (i <= mid) tmp[k++] = v[i++];
	while (j <= r) tmp[k++] = v[j++];
	for (i = l; i <= r; i++) v[i] = tmp[i];
}

void merge(int l, int r)
{
	if (l == r || l > r) return;
	int mid = (l + r) / 2;
	merge(l, mid);
	merge(mid + 1, r);
	if (v[mid] < v[mid + 1]) return;
	merge(l, mid, r);
}

int main()
{
	int n;
	while (scanf("%d", &n) == 1 && n) {
		ans = 0;
		for (int i = 1; i <= n; i++) scanf("%d", &v[i]);
		merge(1, n);
		printf("%lld\n", ans);
	}
}
