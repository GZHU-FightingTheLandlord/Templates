#include <string.h>
#include <queue>
#include <algorithm>
using namespace std;
class Prim{
public:
	void init(int n_ = 0) // 初始化, n为结点数
	{
		cnt = 0;
		n = n_;
		memset(head, -1, (n + 2) * sizeof(int));
	}
	
	void addedge(int u = 0, int v = 0, int w = 0) // 加边 u-v 双向边， 权值为w
	{
		e[++cnt] = edge(v, w, head[u]); head[u] = cnt;
		e[++cnt] = edge(u, w, head[v]); head[v] = cnt;
	}
	
	int solve() // 求最小生成树权值总和
	{
		priority_queue<node> q;
		memset(dis, 0x3f3f3f3f, (n + 2) * sizeof(int));
		
		q.push(node(1, 0));
		dis[1] = 0; ans = 0;
		while (!q.empty())
		{
			int u = q.top().u, tmp = q.top().w;
			q.pop();
			
			if (tmp > dis[u])
				continue;
			ans += dis[u];
			dis[u] = 0;
			for (int i = head[u]; i != -1; i = e[i].next)
			{
				int v = e[i].v, w = e[i].w;
				if (dis[v] > w)
				{
					dis[v] = w;
					q.push(node(v, w));
				}
			}
		}
		return ans;
	}
	
	Prim(){ e = new edge[MAX << 1]; head = new int[MAX]; dis = new int[MAX]; }
private:
	static const int MAX = 1e6 + 5;
	int n, cnt, ans;
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
	};
	edge *e;
	int *head, *dis;
};


