#include <bits/stdc++.h>
using namespace std;

// http://adilet.org/blog/palindromic-tree/

// This one may be better.
struct Palindromic_Tree {
    static const int maxn = 1e5 + 5;

    char s[maxn];
    int n, tot, last, len[maxn], ch[maxn][27], fail[maxn];

    void init() {
        n = tot = last = 0;
        newnode(0); newnode(-1);
        fail[0] = 1; s[n] = -1;
    }

    int newnode(int val) {
        len[tot] = val;
        memset(ch[tot], 0, sizeof ch[tot]);
        return tot++;
    }

    int getfail(int cur) {
        while (s[n] != s[n - len[cur] - 1]) {
            cur = fail[cur];
        }
        return cur;
    }

    void add(char c) {
        c -= 'a'; s[++n] = c;
        int cur = getfail(last);
        if (!ch[cur][c]) {
            int nxt = newnode(len[cur] + 2);
            fail[nxt] = ch[getfail(fail[cur])][c];
            ch[cur][c] = nxt;
        }
        last = ch[cur][c];
    }

    void build(char *c, int length) {
        for (int i = 0; i < length; i++) {
            add(c[i]);
        }
    }
}pt;

// struct Palindromic_Tree {
//     static const int maxn = 2e6 + 5;
//     static const int char_db = 10;

//     int tot, len[maxn], fail[maxn], ch[maxn][char_db];
//     long long sum[maxn];

//     // even root -> 0, odd root -> 1

//     inline int newnode(int val) {
//         sum[tot] = 0, len[tot] = val;
//         memset(ch[tot], 0, sizeof ch[tot]);
//         return tot++;
//     }

//     void init() {
//         tot = 0;
//         newnode(0); newnode(-1);
//         fail[0] = 1, fail[1] = 0;
//     }

//     int getfail(char *s, int cur, int i) {
//         while (i - len[cur] - 1 < 0 || s[i] != s[i - len[cur] - 1]) {
//             cur = fail[cur];
//         }
//         return cur;
//     }

//     void build(char *s, int n) {
//         int cur = 1;
//         for (int i = 0; i < n; i++) {
//             cur = getfail(s, cur, i);
//             if (!ch[cur][s[i] - '0']) {
//                 int nxt = newnode(len[cur] + 2);
//                 fail[nxt] = ch[getfail(s, fail[cur], i)][s[i] - '0'];
//                 ch[cur][s[i] - '0'] = nxt;
//             }
//             cur = ch[cur][s[i] - '0'];
//         }
//     }
// }pt;