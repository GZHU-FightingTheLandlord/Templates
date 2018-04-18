#include <algorithm>
#include <string.h>
using namespace std;
class union_find{ // 并查集引用
public:
	union_find(int _len = 0)
	{
		opn=0;
		len=_len+5;
		data=NULL;
		data=new int[len];
		//if(data==NULL)throw std::runtime_error("ERROR!");
		for(int i=1;i<len;i++)data[i]=i;
	}
	void init()
	{
		for (int i = 0; i < len; i++)
			data[i] = i;
	}
	int find(int x)
	{
		//if(x>len)throw std::runtime_error("ERROR!");
		return x==data[x]?x:data[x]=find(data[x]);
	}
	void merge(int x,int y)
	{
		x=find(x);
		y=find(y);
		if(x!=y)
		{
			data[x]=y;
			opn++;
		}
	}
	int opn;
private:
	int *data;
	int len;
};

class Kruskal{
public:
	void init(int n_ = 0) // 初始化, n_ 为结点数
	{
		n = n_; cnt = 0; ans = 0;
		p.init();
	}
	
	void addedge(int u, int v, int w) // 加边u-v
	{
		e[++cnt] = edge(u, v, w);
	}
	
	int solve() // 求最小生成树权值总和
	{
		sort(e + 1, e + 1 + cnt);
		for (int i = 1; i <= cnt; i++)
		{
			int a = p.find(e[i].u), b = p.find(e[i].v);
			if (a != b)
			{
				ans += e[i].w;
				p.merge(a, b);
			}
		}
		return ans;
	}
	
	Kruskal(){ p = union_find(MAX); e = new edge[MAX];}
private:
	static const int MAX = 1e6 + 5;
	struct edge{
		int u, v, w;
		edge(int u_ = 0, int v_ = 0, int w_ = 0):u(u_), v(v_), w(w_){}
		bool operator<(const edge b)const
		{
			return w < b.w;
		}
	};
	edge *e;
	union_find p;
	int ans, cnt, n;
};
