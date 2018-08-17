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
	// Initial
	void init(int n = -1) {
		(n != -1) ? (cnt = N = n) : (cnt = N);
		for (int i = 0; i <= N; i++) Root[i] = i;
	}
	// Get Ancestor
	int find(int x) {
		return (x == Root[x]) ? x : Root[x] = find(Root[x]);
	}

	bool Union(int a, int b) {
		if (!same(a, b)) {
			cnt--; // Trees in forest
			Root[Root[b]] = Root[a]; // Link b to a
			return true;
		}
		return false;
	}
	// Same Ancestor?
	inline bool same(int a, int b) {
		return find(a) == find(b);
	}
};