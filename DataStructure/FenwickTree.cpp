#include <algorithm>
#include <string.h>
using namespace std;

namespace Fenwick {
#define MAXN 100005
#define lowbit(i)   (i & (-i))
    int _Size;
    int _Sum[MAXN];
    void init(int n = 0) {
        _Size = n;
        memset(_Sum, 0, sizeof _Sum);
    }
    void upd(int i, int delta) {
        for (; i <= _Size; i += lowbit(i)) _Sum[i] += delta;
    }
    int getsum(int i) {
        int ret = 0;
        for (; i > 0; i -= lowbit(i)) ret += _Sum[i];
        return ret;
    }
    int query(int l, int r) {
        return getsum(r) - getsum(l - 1);
    }
#undef MAXN
}
using Fenwick::upd;
using Fenwick::query;
