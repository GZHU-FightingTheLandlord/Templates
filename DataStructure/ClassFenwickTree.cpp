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

    // 询问区间[1, i]的元素总和
    int query(int i)
    {
        int res = 0;
        while (i > 0) {
            res += sum_[i];
            i -= lowbit(i);
        }
        return res;
    }

    FenwickTree(int _n = 0) : sum_(_n + 1, 0), n(_n) {}
private:
    int n;
    vector<int> sum_;

    inline int lowbit(int x)
    {
        return x & (-x);
    }
};
