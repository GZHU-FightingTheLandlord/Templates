#include <algorithm>
#include <string.h>
#include <queue>

struct Aho_Corasick {
    static const int maxn = 5e5 + 5000;
    static const int sigma = 26;

    int tot, son[maxn][sigma], cnt[maxn], fail[maxn];    

    inline void init() {
        tot = cnt[0] = fail[0] = 0;
        memset(son[0], 0, sizeof son[0]);
    }

    inline int trans(int x) {
        return x - 'a';
    }

    void insert(char *s) {
        int now = 0, n = strlen(s);
        for (int i = 0; i < n; i++) {
            int c = trans(s[i]);
            if (!son[now][c]) {
                cnt[++tot] = 0; fail[tot] = 0;
                memset(son[tot], 0, sizeof son[tot]);
                son[now][c] = tot;
            }
            now = son[now][c];
        }
        cnt[now]++;
    }

    std::queue<int> Q;

    void build() {
        while (!Q.empty()) Q.pop();
        for (int i = 0; i < sigma; i++) {
            if (son[0][i]) {
                Q.push(son[0][i]);
            }
        }
        while (!Q.empty()) {
            int u = Q.front(); Q.pop();
            for (int i = 0; i < sigma; i++) {
                if (son[u][i]) {
                    fail[son[u][i]] = son[fail[u]][i];
                    Q.push(son[u][i]);
                }
                else son[u][i] = son[fail[u]][i];
            }
        }
    }

    int solve(char *s) {
        int ret = 0, n = strlen(s);
        for (int i = 0, now = 0, c; i < n; i++, now = son[now][c]) {
            c = trans(s[i]);
            for (int u = son[now][c]; u; u = fail[u]) {
                ret += cnt[u]; cnt[u] = 0;
            }
        }
        return ret;
    }
}Accepted;
