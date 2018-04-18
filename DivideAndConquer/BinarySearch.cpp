#include <bits/stdc++.h>
using namespace std;

int search(int v[], int l, int r, int k)
{
    if (l > r || l < 0 || r >= 100)
        return 0;
    int mid = (l + r) >> 1;
    if (v[mid] < k)
        return search(v, mid + 1, r, k);
    if (v[mid] > k)
        return search(v, l, mid - 1, k);
    return v[mid] == k;
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
