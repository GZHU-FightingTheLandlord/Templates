#include <string.h>

struct Trie {
    static const int maxn = 5e5 + 5;
    static const int sigma = 26;
    int tot, ch[maxn][sigma], cnt[maxn];
    bool isend[maxn];

    inline int newnode() {
        ++tot;
        memset(ch[tot], 0, sizeof ch[tot]);
        cnt[tot] = 0, isend[tot] = false;
        return tot;
    }

    inline void init() { tot = -1, newnode(); }

    inline int trans(char c) {
        // edit
        // Example: return c - 'a';
    }

    void insert(char *str, int len) {
        int now = 0;
        for (int i = 0; i < len; i++) {
            int c = trans(str[i]);
            if (!ch[now][c]) {
                ch[now][c] = newnode();
            }
            cnt[now]++;
            now = ch[now][c];
        }
        cnt[now]++;
        isend[now] = true;
    }

    bool find(char *str, int len) {
        int now = 0;
        for (int i = 0; i < len; i++) {
            if (!ch[now][c]) return false;
            now = ch[now][c];
        }
        return isend[now];
    }
}trie;