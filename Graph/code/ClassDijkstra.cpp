#include <algorithm>
#include <string.h>
#include <vector>
#include <queue>
#include <stack>
using namespace std;

class Dijkstra{
public:
	void init(int n_ = 0) // 初始化, n_为结点数
	{
		cnt = 0;
		n = n_;
		memset(head, -1, sizeof head);
	}
	
	void addedge(int u, int v, int w) // 单向加边u->v, w为权值
	{
		e[cnt] = edge(v, w, head[u]); head[u] = cnt++;
	}
	
	int solve(int s, int t) // 求s->t最短路程
	{
		priority_queue<node> q;
		memset(dis, 0x3f3f3f3f, sizeof dis);
		memset(vis, 0, sizeof vis);
        memset(pre, -1, sizeof pre);
		
		q.push(node(s, 0));
		dis[s] = 0;
		while (!q.empty())
		{
			int u = q.top().u;
			q.pop();
			if (u == t)
				return dis[t];
			if (vis[u] != 0)
				continue;
			vis[u] = 1;
			for (int i = head[u]; i != -1; i = e[i].next)
			{
				int v = e[i].v, w = e[i].w;
				if (!vis[v] && dis[v] > dis[u] + w)
				{
					dis[v] = dis[u] + w;
                    pre[v] = u;
					q.push(node(v, dis[v]));
				}
			}
		}
		return -1;
	}

    vector<int> getpath(int s, int t) // 求s->t最短路径
    {
        solve(s, t);
        stack<int> st;
        vector<int> ans;
        int i = t;
        while (i != -1)
        {
            st.push(i);
            i = pre[i];
        }
        while (!st.empty())
        {
            ans.push_back(st.top());
            st.pop();
        }
        return ans;
    }
    
	Dijkstra(){}
private:
	static const int MAX = 1e5 + 5;
	int cnt, n;
	struct node{
		int u, w;
		node(int u_ = 0, int w_ = 0):u(u_), w(w_){}
		bool operator<(const node b)const
		{
			return w > b.w;
		}
	};
	struct edge{
		int v, w;
		int next;
		edge(int v_ = 0, int w_ = 0, int n_ = -1):v(v_), w(w_), next(n_){}
	}e[MAX];
	int dis[MAX];
	int head[MAX];
    int pre[MAX];
	bool vis[MAX];
};


