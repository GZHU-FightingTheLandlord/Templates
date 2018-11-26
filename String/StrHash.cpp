#include <vector>
using namespace std;
typedef unsigned long long ull;
struct strhash {
    static const ull seed = 1313131;
    vector<ull> h, p;
    strhash(int n = 0) : h(n + 5, 0), p(n + 5, 0) {}
    void init(char *s) {
        for (int i = 0; s[i]; i++) {
            if (i) h[i] = h[i - 1] * seed, p[i] = p[i - 1] * seed;
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