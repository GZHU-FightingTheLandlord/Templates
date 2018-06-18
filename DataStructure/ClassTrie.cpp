#include <algorithm>
#include <string.h>
#include <vector>
using std::vector;

/*
usage:

    Define: Trie T(MaxSize);

    Initial: T.init();

    Push in: T.push(Str, StrLen);

    Check: T.check(Str, StrLen);
*/

struct Trie {

    struct Node {
        bool is_end;
        int next[30];
        Node() { is_end = 0; memset(next, 0, sizeof next); }
    };

    int cnt;
    vector<Node> Vn;

    // Max Size!!!!!
    Trie(int Size = 0) : Vn(Size + 5) {
        cnt = 0;
    }

    void init()
    {
        cnt = 0;
        Vn[0] = Node();
    }

    void push(const char *arr, const int& len)
    {
        int now = 0;
        for (int i = 0; i < len; i++) {
            if (Vn[now].next[arr[i] - 'a'] == 0) {
                Vn[++cnt] = Node();
                Vn[now].next[arr[i] - 'a'] = cnt;
            }
            now = Vn[now].next[arr[i] - 'a'];
        }
        Vn[now].is_end = true;
    }

    bool check(const char *arr, const int& len)
    {
        int now = 0;
        for (int i = 0; i < len; i++) {
        if (Vn[now].next[arr[i] - 'a'] == 0) {
                return false;
            }
            now = Vn[now].next[arr[i] - 'a'];
        }
        return Vn[now].is_end;
    }
};
