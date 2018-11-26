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
        int u = 0;
        for (int i = 0; i < len; i++) {
            int c = trans(str[i]);
            if (!ch[u][c]) ch[u][c] = newnode();
            cnt[u]++;
            u = ch[u][c];
        }
        cnt[u]++;
        isend[u] = true;
    }

    bool find(char *str, int len) {
        int u = 0;
        for (int i = 0; i < len; i++) {
            char c = trans(str[i]);
            if (!ch[u][c]) return false;
            u = ch[u][c];
        }
        return isend[u];
    }
}trie;

struct ZeroOneTrie {
    static const int N = 1e5 + 5; // total number of integer
    static const int maxn = N * 31; // 31 bit for int

    int tot;
    int ch[maxn][2], cnt[maxn];

    void init() { tot = ch[0][0] = ch[0][1] = cnt[0] = 0; }
    
    // insert integer x into Trie
    void insert(int x) {
        int u = 0; cnt[u]++;
        for (int i = 30; ~i; --i) {
            int v = (x >> i) & 1;
            if (!ch[u][v]) {
                int tmp = ++tot;
                ch[tmp][0] = ch[tmp][1] = cnt[tmp] = 0;
                ch[u][v] = tmp;
            }
            u = ch[u][v]; ++cnt[u];
        }
    }

    // erase integer x from Trie(if x does not exist, it will cause some ERROR)
    void erase(int x) {
        int u = 0; cnt[u]--;
        for (int i = 30; ~i; --i) {
            int v = (x >> i) & 1;
            u = ch[u][v]; --cnt[u];
        }
    }

    // max(x xor a_i), 1 \leq i \leq a_size
    int MaxXor(int x) {
        int ret = 0, u = 0;
        for (int i = 30; ~i; --i) {
            int v = (x >> i) & 1;
            if (ch[u][v ^ 1] && cnt[ch[u][v ^ 1]] > 0) {
                ret |= (v ^ 1) << i;
                u = ch[u][v ^ 1];
            }
            else {
                ret |= v << i;
                u = ch[u][v];
            }
        }
        return ret ^ x;
    }
};