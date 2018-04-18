#include <string.h>
#include <vector>
#include <queue>
#include <stack>
using namespace std;

class spfa{
public:
	void init(int n_ = 0) // 初始化， n为结点数
	{
		cnt = 0;
		n = n_;
		memset(head, -1, sizeof head);
	}
	
	void addedge(int u, int v, int w) // 单向加边,u -> v， 权值为w
	{
		e[cnt] = edge(v, w, head[u]); head[u] = cnt++;
	}
	
	int solve(int s, int t) // 求s -> t的最短路程
	{
		queue<int> q;
		memset(inq, 0, sizeof inq);
		memset(dis, 0x3f3f3f3f, sizeof dis);
        memset(pre, -1, sizeof pre);
		
		dis[s] = 0; inq[s] = 1;
		q.push(s);
		while (!q.empty())
		{
			int u = q.front();
			q.pop(); inq[u] = 0;
			
			for (int i = head[u]; i != -1; i = e[i].next)
			{
				int v = e[i].v, w = e[i].w;
				if (dis[v] > dis[u] + w)
				{
					dis[v] = dis[u] + w;
                    pre[v] = u;
					if (!inq[v])
						q.push(v), inq[v] = 1;
				}
			}
		}
		return dis[t] == 0x3f3f3f3f ? -1 : dis[t];
	}

    vector<int> getpath(int s, int t) // 求s -> t的最短路径
    {
        vector<int> p;
        solve(s, t);
        stack<int> ss;
        int i = t;
        while (i != -1)
        {
            ss.push(i);
            i = pre[i];
        }
        while (!ss.empty())
        {
            p.push_back(ss.top());
            ss.pop();
        }
        return p;
    }
    
	spfa(){}
private:
	static const int MAX = 1e5 + 5;
	int cnt, n;
	struct edge{
		int v, w;
		int next;
		edge(int v_ = 0, int w_ = 0, int n_ = -1):v(v_), w(w_), next(n_){}
	}e[MAX];
	int dis[MAX];
	int head[MAX];
    int pre[MAX];
	bool inq[MAX];
};
