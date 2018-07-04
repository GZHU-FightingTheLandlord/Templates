#include <algorithm>
using namespace std;

int v[1005];
void Qsort(int l, int r)
{
    if (l >= r) return;

    swap(v[l], v[(rand() % (r - l + 1)) + l]);

    int i = l, j = r + 1;
    while (1) {
        while (v[l] > v[++i]);
        while (v[l] < v[--j]);
        if (i >= j) break;
        swap(v[i], v[j]);
    }
    swap(v[l], v[j]);

    Qsort(l, j - 1);
    Qsort(j + 1, r);
}