#include <algorithm>
#include <string.h>
using namespace std;

const int maxn = 110005;

int p[maxn << 1];
char str[maxn << 1];

int manacher(char *s, int n) {
    str[0] = '$'; str[1] = '#';

    for (int i = 0; i < n; i++) {
        str[(i << 1) + 2] = s[i];
        str[(i << 1) + 3] = '#';
    }
    n = (n + 1) << 1;
    str[n] = 0;

    int ret = 0, mx = 0, pos;
    for (int i = 1; i < n; i++) {
        p[i] = mx > i ? min(p[(pos << 1) - i], mx - i) : 1;
        
        while (str[i - p[i]] == str[i + p[i]]) p[i]++;
        
        if (p[i] + i > mx) mx = p[i] + i, pos = i;

        ret = max(ret, p[i]);
    }
    return ret - 1;
}

// ****************************************************

// index start from one
int solve(char *s, int n) {
    int mx = 1; s[0] = '$'; s[++n] = 0;
    for (int i = 0, p, q; i < n; i++) {
        for (q = i; s[i + 1] == s[i]; ++i);
        for (p = i; s[q - 1] == s[p + 1]; --q, ++p);
        mx = max(mx, p - q + 1);
    }
    return mx;
}