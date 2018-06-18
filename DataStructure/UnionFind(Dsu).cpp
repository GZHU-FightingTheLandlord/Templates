#include <algorithm>
#include <vector>
using std::vector;

// Disjoint-set union

struct Dsu {
	int N, cnt;
	vector<int> Root;

	// N -> MaxNum
	Dsu(int num = 0) : N(num), Root(num + 5), cnt(num) {
		for (int i = 0; i < num; i++) Root[i] = i;
	}

	void init(int n = -1) // Initial
	{
		if (n != -1) N = n;
		cnt = N;
		for (int i = 0; i <= N; i++) Root[i] = i;
	}

	int find(int x) // Get Ancestor
	{
		if (x == Root[x]) return x;
		else return Root[x] = find(Root[x]);
	}

	bool Union(int a, int b)
	{
		if (!same(a, b)) {
			cnt--; // Trees in forest
			Root[Root[b]] = Root[a]; // Link b to a
			return true;
		}
		return false;
	}

	bool same(int a, int b)
	{
		a = find(a);
		b = find(b);
		return a == b; // Same Ancestor?
	}
};