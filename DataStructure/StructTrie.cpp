/*
usage:

    Define: Trie T(MaxSize);

    Initial: T.init();

    Push in: T.push(Str, StrLen);

    Check: T.check(Str, StrLen);
*/

#include <algorithm>
#include <string.h>
#include <vector>
using namespace std;

struct Trie {

    struct Node {
        bool ed;
        int next[30];
        Node() { ed = false; memset(next, 0, sizeof next); }
    };

    int cnt;
    vector<Node> v;

    // Max Size!!!!!
    Trie(int Size = 0) : v(Size + 5) { cnt = 0; }

    void init()
    {
        cnt = 0;
        v[0] = Node();
    }

    void push(const char *arr, const int& len)
    {
        int now = 0;
        for (int i = 0; i < len; i++) {
            if (v[now].next[arr[i] - 'a'] == 0) {
                v[++cnt] = Node();
                v[now].next[arr[i] - 'a'] = cnt;
            }
            now = v[now].next[arr[i] - 'a'];
        }
        v[now].ed = true;
    }

    bool check(const char *arr, const int& len)
    {
        int now = 0;
        for (int i = 0; i < len; i++) {
        if (v[now].next[arr[i] - 'a'] == 0) {
                return false;
            }
            now = v[now].next[arr[i] - 'a'];
        }
        return v[now].ed;
    }
};
