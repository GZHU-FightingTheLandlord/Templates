#include <algorithm>
#include <vector>
using namespace std;

class FenwickTree{
public:
    // i处元素增加delta
    void update(int i, int delta)
    {
        while (i <= n) {
            sum_[i] += delta;
            i += lowbit(i);
        }
    }

    // 询问区间[l, r]的元素总和
    int query(int l, int r)
    {
        return getres(r) - getres(l - 1);
    }

    FenwickTree(int _n = 0) : sum_(_n + 1, 0), n(_n) {}
private:
    int n;
    vector<int> sum_;

    int getres(int i)
    {
        int res = 0;
        while (i > 0) {
            res += sum_[i];
            i -= lowbit(i);
        }
        return res;
    }

    inline int lowbit(int x)
    {
        return x & (-x);
    }
};
