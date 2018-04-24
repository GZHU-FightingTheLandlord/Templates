/*
    待优化
*/
#include <bits/stdc++.h>
using namespace std;
/* lower_bound近似    */
int search(int v[], int l, int r, int k)
{
	int mid;
	while (l < r)
	{
		mid = (l + r) >> 1;
		if (v[mid] < k)
			l = mid + 1;
		else
			r = mid;
	}
	return l;
}

int main()
{
    int n;
    int v[100];
    scanf("%d", &n);
    for (int i = 0; i < n; i++)
        scanf("%d", &v[i]);
    sort(v, v + n);
    int q;
    scanf("%d", &q);
    while (q--)
    {
        int k;
        scanf("%d", &k);
        int tmp = search(v, 0, n - 1, k);
        printf("%d\n", tmp);
        /*
            int tmp = lower_bound(v, v + n, k) - v;
        */
    }
    return 0;
}
