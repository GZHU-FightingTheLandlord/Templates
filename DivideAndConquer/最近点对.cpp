#include <iostream>
#include <algorithm>
#include <math.h>
using namespace std;

struct point {
    double x, y;
    point() {};
}v[100005], tmp[100005];

inline double sqr(const double &x) { return x * x; }

inline double dist(const point &a, const point &b)
{
    return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
}

double Div(int l, int r)
{
    if (l >= r) {
        return 1e308;
    }
    if (l + 1 == r) {
        return dist(v[l], v[r]);
    }
    int mid = (l + r) >> 1;
    double ret = min(Div(l, mid), Div(mid + 1, r));

    int k = 0;
    for (int i = mid; i >= l && !(fabs(v[i].x - v[mid].x) > ret); i--) {
        tmp[k++] = v[i];
    }
    for (int i = mid + 1; i <= r && !(fabs(v[i].x - v[mid].x) > ret); i++) {
        tmp[k++] = v[i];
    }

    sort(tmp, tmp + k, [](const point& a, const point& b) {
        return a.y < b.y;
    });

    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            if (fabs(tmp[i].y - tmp[j].y) > ret) {
                break;
            }
            ret = min(ret, dist(tmp[i], tmp[j]));
        }
    }
    return ret;
}

double solve(int st, int n)
{
    sort(v + st, v + st + n);
    return Div(st, st + n - 1);
}

int main()
{
    ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);

    int n;
    while (cin >> n) {
        for (int i = 0; i < n; i++) {
            cin >> v[i].x >> v[i].y;
        }
        cout << solve(0, n) << endl;
    }
}