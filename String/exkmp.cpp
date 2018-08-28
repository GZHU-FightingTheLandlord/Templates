// nxt[i]: t[i...m-1]与t[0...m-1]的最长公共前缀
// extend[i]: s[i...n-1]与t[0...m-1]的最长公共前缀

const int maxn = 1e6 + 5;
int nxt[maxn], extend[maxn];
void exkmp(char *s, int n, char *t, int m) {
    int j = 0, k = 1;
    for (; j + 1 < m && t[j] == t[j + 1]; ++j);
    nxt[0] = m, nxt[1] = j;

    for (int i = 2; i < m; i++) {
        int p = nxt[k] + k - 1, L = nxt[i - k];
        if (p + 1 - i - L > 0) {
            nxt[i] = L;
        }
        else {
            for (j = (p - i + 1 > 0) ? (p - i + 1) : 0; i + j < m && t[i + j] == t[j]; ++j);
            // blank
            nxt[i] = j, k = i;
        }
    }

    j = k = 0;
    for (; j < n && j < m && t[j] == s[j]; ++j);
    extend[0] = j;

    for (int i = 1; i < n; i++) {
        int p = extend[k] + k - 1, L = nxt[i - k];
        if (p + 1 - i - L > 0) {
            extend[i] = L;
        }
        else {
            for (j = (p - i + 1 > 0) ? (p - i + 1) : 0; i + j < n && j < m && s[i + j] == t[j]; ++j);
            // blank
            extend[i] = j, k = i;
        }
    }
}