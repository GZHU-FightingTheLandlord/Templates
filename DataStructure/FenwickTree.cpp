#include <stdio.h>
#include <algorithm>
#include <vector>
using std::vector;

template <typename T> struct Fenwick{
	int n;
	vector<T> sum_;
    
	#define lowbit(x) (x & (-x))
		
	Fenwick(int n_ = 0) : n(n_), sum_(n + 5, 0) {}
    
	void upd(int i, T delta)
	{
        for (; i <= n; i += lowbit(i)) sum_[i] += delta;
	}
	T getsum(int i)
    {
		T res = 0;
        for (; i > 0; i -= lowbit(i)) res += sum_[i];
		return res;
	}
    
	T query(int l, int r)
    {
		return getsum(r) - getsum(l - 1);
    }
};
