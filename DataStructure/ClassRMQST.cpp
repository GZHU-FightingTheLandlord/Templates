#include <string.h>
#include <algorithm>
#include <math.h>
using namespace std;

class RMQST{

public:

    void init(int n_, int v[])
    {
        n = n_;
        memcpy(Min[0], v, sizeof Min[0]);
        solve();
    }

    int query(int l, int r) // O(1) 询问
    {
        int mid = (int)((log(r - l + 1)) / (log(2.0)));
        /* 按需选择返回方式     */
        return max(Max[l][mid], Max[r - (1 << mid) + 1][mid]);
        //return min(Min[l][mid], Min[r - (1 << mid) + 1][mid]);    
    }

    RMQST(){ n = 0; memset(Min, 0, sizeof Min); memset(Max, 0, sizeof Max); }

private:

    static const int MAX = 1e5 + 5;

    int n;
    int Min[MAX][25];
    int Max[MAX][25];

    void solve()  // O(NlogN) 建表
    {
        int l = (int)((log(n)) / (log(2.0)));
        for (int j = 1; j <= l; j++)
        {
            for (int i = 1; i + (1 << j) - 1 <= n; i++)
            {
                Min[i][j] = min(Min[i][j - 1], Min[i + (1 << (j - 1))][j - 1]);
                Max[i][j] = max(Max[i][j - 1], Max[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

};

