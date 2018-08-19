#include <algorithm>
#include <cmath>
using namespace std;
const int MAX = 10005;
const double eps = 1e-8;

/*
    O(nlogn)
    p->原点集, s->凸包
    Convex::Graham(n); // 返回值为凸包中点数
*/

namespace Convex {
    struct Point {
        double x, y;
        Point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    }p[MAX], s[MAX];
    template <typename T> inline T sqr(T x) { return x * x; }
    inline int dcmp(double x) {
        if (fabs(x) <= eps) return 0;
        return x < 0 ? -1 : 1;
    }
    inline double dist(Point a, Point b) {
        return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
    }
    inline double cross(Point a, Point b, Point c) {
        return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
    }
    bool cmp(const Point& a, const Point& b) {
        int c = dcmp(cross(p[0], a, b));
        if (c < 0) return false;
        if (c > 0) return true;
        return dist(p[0], a) < dist(p[0], b);
    }
    int Graham(int n) {
        for (int i = 0; i < n; i++) {
            if (dcmp(p[i].y - p[0].y) < 0) {
                swap(p[i], p[0]);
            }
            else if (!dcmp(p[i].y - p[0].y) && dcmp(p[i].x - p[0].x) < 0) {
                swap(p[i], p[0]);
            }
        }
        sort(p + 1, p + n, cmp);
        s[0] = p[0], s[1] = p[1];
        int tot = 1;
        for (int i = 2; i < n; i++) {
            while (tot > 0 && dcmp(cross(s[tot - 1], s[tot], p[i])) <= 0) {
                --tot;
            }
            s[++tot] = p[i];
        }
        return tot + 1;
    }
}