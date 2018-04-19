#include <bits/stdc++.h>
using namespace std;

const int MAX = 505;

int n, m;
vector<int> e[MAX];
int deg[MAX];
vector<int> ans;

/*   dfs实现， 判环用， 需对所有为vis点调用一次

int vis[MAX];
bool tps(const int& u)
{
    vis[u] = 1;
    for (auto v : e[u])
        if (vis[v] == 1 || (!vis[v] && !tps(v)))
            return 0;
    ans.push_back(u);
    vis[u] = 2;
    return 1;
}

*/

bool tps()
{
    /* 优先队列不是必须的    */
    priority_queue<int, vector<int>, greater<int> > q;

    /* 所有入度为0的点进队   */
    for (int i = 1; i <= n; i++)
        if (!deg[i])
            q.push(i);

    while (!q.empty())
    {
        int u = q.top();
        q.pop();

        ans.push_back(u);
        for (auto v : e[u])
        {
            deg[v]--;
            if (!deg[v]) // 入度为0的点入队
                q.push(v);
        }
    }

    /* 若有入度不为0的点， 则图中存在环    */
    for (int i = 1; i <= n; i++)
        if (deg[i])
            return 0;
    return 1;
}

int main()
{
    scanf("%d%d", &n, &m);
    while (m--)
    {
        int u, v;
        scanf("%d%d", &u, &v);
        e[u].push_back(v);  // 建边
        deg[v]++;           // v入度加1
    }
    if (!tps())
        printf("No Answer\n");
    else
    {
        for (auto i : ans)
            printf("%d ", i);
        printf("\n");
    }
    return 0;
}
