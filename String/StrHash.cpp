#include <vector>
using namespace std;
typedef unsigned long long ull;
struct strhash {
    vector<ull> h, p;
    strhash(int n = 0) : h(n + 5, 0), p(n + 5, 0) {}
    void init(char *s) {
        for (int i = 0; s[i]; i++) {
            if (i) h[i] = h[i - 1] * 131, p[i] = p[i - 1] * 131;
            else p[i] = 1;
            h[i] += s[i] - 'a' + 1;
        }
    }
    inline ull gethash(int l, int r) {
        ull ret = h[r];
        if (l) ret -= h[l - 1] * p[r - l + 1];
        return ret;
    }
};