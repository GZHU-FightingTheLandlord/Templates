#include <algorithm>
#include <vector>
#include <queue>
using namespace std;

struct TSDFS {
#define szz(s) (int)s.size()
#define all(a) a.begin(),a.end()
    int N;
    bool circle;
    vector<vector<int>> e;
    vector<int> seq, vis;

    TSDFS(int n = 0) : N(n), e(n + 5), seq(n + 5), vis(n + 5) { init(); }

    void init(int n = -1)
    {
        if (n != -1) N = n;
        for (int i = 0; i <= N; i++) {
            e[i] = vector<int>();
            vis[i] = 0;
        }
    }

    void addedge(int u, int v)
    {
        e[u].push_back(v);
    }


    void dfs(int u)
    {
        vis[u] = 1;
        for (int i = 0; i < szz(e[u]); i++) {
            int& v = e[u][i];
            (vis[v] == 1 && (circle = true)) || (vis[v] == 0 && (dfs(v), 1));
        }
        vis[u] = 2;
        seq.push_back(u);
    }


    void solve()
    {
        circle = false; seq.clear();
        for (int i = 1; i <= N; i++) {
            if (!vis[i]) {
                dfs(i);
            }
        }
        reverse(all(seq));
    }
};

struct TSBFS {
#define szz(a) (int)a.size()
#define all(a) a.begin(),a.end()
    int N;
    bool circle;
    vector<vector<int>> e;
    vector<int> seq, deg;

    TSBFS(int n = 0) : N(n), circle(false), e(n + 5), seq(n + 5), deg(n + 5) {}

    void init(int n = -1) 
    {
        if (n != -1) N = n;
        for (int i = 0; i <= N; i++) {
            e[i] = vector<int>();
            deg[i] = 0;
        }
    }

    void addedge(int u, int v)
    {
        e[u].push_back(v);
        deg[v]++;
    }

    void solve()
    {
        queue<int> Q;
        for (int i = 1; i <= N; i++) {
            if (deg[i] == 0) {
                Q.push(i);
            }
        }

        seq.clear(); circle = false;
        while (!Q.empty()) {
            int u = Q.front(); Q.pop();

            seq.push_back(u);

            for (int i = 0; i < szz(e[u]); i++) {
                int& v = e[u][i];
                if ((--deg[v]) == 0) {
                    Q.push(v);
                }
            }
        }

        for (int i = 1; i <= N; i++) {
            circle = circle || (deg[i] > 0);
        }
    }
};