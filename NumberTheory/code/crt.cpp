/*
  链接：http://poj.org/problem?id=2115
*/

#include<iostream>
using namespace std;
typedef long long ll;

ll mod(ll x, ll p)
{
    return ((x % p) + p) % p;
}

ll inv, y;
ll exgcd(long a, long b)
{
    if (b == 0)
    {
        inv = 1;
        y = 0;
        return a;
    }
    else
    {
        long ret = exgcd(b, a % b);
        long tmp = inv;
        inv = y;
        y = tmp - a / b * y;
        return ret;
    }
}

signed main()
{
#ifndef ONLINE_JUDGE
    freopen("test.txt", "r", stdin);
#endif
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);
    //TODO:
    int T, A, R;
    while (cin >> T)
    {
        bool flag = false;
        ll a = 1, r = 0;
        while (T--)
        {
            cin >> A >> R;
            const ll t = exgcd(a, A);
            // t为a和A的gcd
            if ((R - r) % t != 0)
            {
                flag = true;
            }
            // inv为a和A的inv
            r += mod(inv * (R - r) / t, A / t) * a;
            a *= A / t;
        }
        if (flag)
        {
            cout << -1 << endl;
        }
        else
        {
            cout << mod(r, a) << endl;
        }
    }
}
