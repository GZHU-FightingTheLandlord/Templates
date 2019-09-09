### min_25用途
1. 筛出质数
2. 求出所有的$\sum\limits_{i=2}^{\lfloor\frac{n}{x}\rfloor}[i\text{是质数}]f(i)$
3. 求积性函数前缀和

### 前提
1. 当$i$为质数时$f(i)$需要是一个多项式。
2. 对于求积性函数前缀和而言$f(p^k)$需要快速求出，求多个值的时候一般不适用min_25筛。

### 计算

#### 处理质数

$$g(a,b)=\sum_{i=2}^a [\text{$i$ 是质数 或 $pmin_i>prime_b$}] * i^k$$

需要求每一个

$$g(\lfloor\frac{n}{x}\rfloor,\infty)=\sum\limits_{i=2}^{\lfloor\frac{n}{x}\rfloor}[i\text{是质数}]f(i)$$

那么有

$$g(a,b)=
\begin{cases}
g(a,b-1), & a<prime_b^2\\
g(a,b-1)-prime_b^k\left(g(\lfloor \frac{a}{prime_b} \rfloor,b-1)-g(prime_{b-1},b-1)\right), & a\ge prime_b^2
\end{cases}$$

滚动数组叠上去即可求出。

#### 计算前缀和

$$S(a,b)=\sum_{i=2}^a [pmin_i\ge prime_b]f(i)$$

前缀和即

$$\sum\limits_{i=1}^{n}f(i)=S(n,1)+f(1)$$

那么有

$$S(a,b)=
\begin{cases}
0, & a<prime_b\\
\\
g(a,\infty)-g(prime_{b-1},\infty)+ \\
\sum\limits_{i=b}^{\infty} \sum\limits_{t\ge 1, prime_i^{t+1}\le a}\left( S(\lfloor \frac{a}{prime_i^t} \rfloor,i+1) * f(prime_i^t) + f(prime_i^{t+1}) \right), &a\ge prime_b
\end{cases}$$

无需记忆化，递归求解。

### 防忘代码

题目链接：[欧拉函数之和](https://www.51nod.com/Challenge/Problem.html#problemId=1239)

```cpp
#pragma GCC optimize("-O3")
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
void Main();
#ifdef ConanYu
#include "local.hpp"
#else
#define debug(...) do { } while(0)
int main() {
  ios::sync_with_stdio(false), cin.tie(0), cout.tie(0);
  Main();
}
#endif

const int MOD = 1e9 + 7;
int fpow(int a, int b) {
  int ans = 1;
  for(; b; b >>= 1, a = 1ll * a * a % MOD) {
    if(b & 1) ans = 1ll * ans * a % MOD;
  }
  return ans;
}

const int N = 1e5 + 10, INV2 = fpow(2, MOD - 2);
ll n, a[N << 1];
int sn, cnt, id, prime[N], g[N << 1], h[N << 1];
int idx(ll x) {
  return x <= sn ? x : id - n / x + 1;
}

void sub(int &a, int b) {
  a -= b;
  if(a < 0) a += MOD;
}

int solve(long long a, int b) {
  if(a < prime[b]) return 0;
  int ans = g[idx(a)];
  sub(ans, g[idx(prime[b - 1])]);
  for(int i = b; i <= cnt && a / prime[i] >= prime[i]; i++) {
    ll k = prime[i], m = k - 1;
    while(a / k >= prime[i]) {
      ans += 1ll * solve(a / k, i + 1) * m % MOD;
      if(ans >= MOD) ans -= MOD;
      m = m * prime[i] % MOD;
      ans += m;
      if(ans >= MOD) ans -= MOD;
      k *= prime[i];
    }
  }
  return ans;
}
void Main() {
  cin >> n;
  sn = sqrt(n);
  cnt = id = 0;
  for(ll i = 1; i <= n; i = a[id] + 1) {
    a[++id] = n / (n / i);
    g[id] = a[id] % MOD * ((a[id] + 1) % MOD) % MOD * INV2 % MOD;
    sub(g[id], 1);
    h[id] = (a[id] - 1) % MOD;
  }
  for(int i = 2; i <= sn; i++) {
    if(h[i] != h[i - 1]) {
      prime[++cnt] = i;
      for(int j = id; a[j] / i >= i; j--) {
        const int ta = idx(a[j] / i), tb = idx(prime[cnt - 1]);
        int A = g[ta], B = h[ta];
        sub(A, g[tb]), sub(B, h[tb]);
        sub(g[j], 1ll * i * A % MOD);
        sub(h[j], B);
      }
    }
  }
  for(int i = 1; i <= id; i++) {
    g[i] -= h[i];
    if(g[i] < 0) g[i] += MOD;
  }
  cout << ((solve(n, 1) + 1) % MOD) << "\n";
}
```

### 参考来源
[Min-25筛学习笔记](https://cmxrynp.github.io/2018/12/03/Min-25%E7%AD%9B%E5%AD%A6%E4%B9%A0%E7%AC%94%E8%AE%B0/)