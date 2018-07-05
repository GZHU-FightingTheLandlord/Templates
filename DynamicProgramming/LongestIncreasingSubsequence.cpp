// O(nlogn)

#include <iostream>
#include <algorithm>
using namespace std;
   
const int MAX = 500005;
typedef long long ll;

ll v[MAX];
ll stack[MAX];
   
int upb(int l, int r, int k)
{
    int mid;
    while (l < r) {
        mid = (l + r) >> 1;
        if (stack[mid] > k) r = mid;
        else l = mid + 1;
    }
    return r;
}

void solve()
{
    int n; cin >> n;
    for (int i = 1; i <= n; i++) cin >> v[i];

    int top = 0;
    for (int i = 1; i <= n; i++) {
        if (top == 0) stack[++top] = v[i];
        else if (stack[top] <= v[i]) stack[++top] = v[i];
        else {
            int pos = upb(1, top, v[i]);
            stack[pos] = v[i];
        }
    }

    cout << (n - top) << endl;
}
   
int main()
{
    ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
    int times; cin >> times;
    while (times--) solve();
}