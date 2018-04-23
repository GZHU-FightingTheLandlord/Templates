#include <algorithm>
#include <vector>
#include <string.h>
using namespace std;

class Trie{
public:
	void init(int k_ = 1) // 初始化， k_为预估需求空间
	{
		cnt = 0;
		v.resize(k_);
		v[0] = node();
	}
	
	void put(char t[], int len) // 把t[]放入Trie， len为t[]长度
	{
		int now = 0;
		for (int i = 0; i < len; i++)
		{
			if (v[now].next[t[i] - 'a'] == 0)
			{
				v[now].next[t[i] - 'a'] = ++cnt;
				if (cnt == v.size())
					v.push_back(node());
				else
					v[cnt] = node();
			}
			now = v[now].next[t[i] - 'a'];
		}
		v[now].ex = 1;
		v[now].cnt++;
	}
	
	int check(char t[], int len) // 查找t[], 找到则返回以其为前缀的单词数量， 否则返回-1
	{
		int now = 0;
		for (int i = 0; i < len; i++)
		{
			if (v[now].next[t[i] - 'a'] == 0)
				return -1;
			now = v[now].next[t[i] - 'a'];
		}
		return v[now].ex == 1 ? v[now].cnt : -1;
	}

	Trie(){}
	
private:
	struct node{
		int ex, cnt;
		int next[30];
		node() { ex = 0; cnt = 0; memset(next, 0, sizeof next);	}
	};
	
	int cnt;
	vector<node> v;
};

