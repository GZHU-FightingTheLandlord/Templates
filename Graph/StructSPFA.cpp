#include <algorithm>
#include <vector>
#include <queue>
using namespace std;

struct SPFA {
    int N;
    vector<vector<pair<int, int>>> e;
    vector<int> dis, inq, pre, cnt;

    SPFA(int n = 0) : N(n), e(n + 5), dis(n + 5), inq(n + 5), pre(n + 5), cnt(n + 5) { init(); }

    void init(int n = -1) 
    {
        if (n != -1) N = n;
        for (int i = 0; i <= N; i++) {
            e[i] = vector<pair<int, int>>();
            dis[i] = 0x3f3f3f3f; inq[i] = 0;
            pre[i] = -1; cnt[i] = 0;
        }
    }

    void addedge(int u, int v, int w)
    {
        e[u].push_back(make_pair(v, w));
        // e[v].push_back(make_pair(u, w));
    }

    bool solve(int st = 1, int ed = -1)
    {
        queue<int> Q;
        Q.push(st); dis[st] = 0; inq[st] = 1;

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop(); inq[u] = 0;

            if (cnt[u] >= N) return false;

            for (int i = 0; i < e[u].size(); i++) {
                int& v = e[u][i].first, w = e[u][i].second;
                if (dis[v] > dis[u] + w) {
                    dis[v] = dis[u] + w;
                    if (!inq[v]) {
                        Q.push(v);
                        pre[v] = u; 
                        cnt[v]++; inq[v] = 1;
                    }
                }
            }
        }
        return ed == -1 || (dis[ed] == 0x3f3f3f3f);
    }
};