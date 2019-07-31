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

题目链接：[P4213 【模板】杜教筛（Sum）](https://www.luogu.org/problem/P4213)

```cpp
#include<cstdio>
#include<math.h>

using namespace std;
#define ll long long

const int N = 46345;
ll g[N << 1];
int T, id, sum[N], h[N << 1];
unsigned cnt, sn, n, id1[N], id2[N], prime[N], a[N << 1];
bool p[N];
inline unsigned Id(unsigned x) {
  return x <= sn ? x : id - n / x + 1;
}
ll SolvePhi(unsigned a, int b) {
  if(a < prime[b]) {
    return 0;
  }
  ll ans = g[Id(a)] - (sum[b - 1] - (b - 1));
  for(unsigned i = b; i <= cnt && prime[i]*prime[i] <= a; ++i) {
    // 这里是强行展开了一层，可能会快一点，因为条件必然满足，事实上可以和下面的写在一起
    ans += (SolvePhi(a / prime[i], i + 1) + prime[i]) * (prime[i] - 1);
    for(unsigned x = prime[i] * prime[i], f = x - prime[i]; (ll)x * prime[i] <= a; x = x * prime[i], f *= prime[i]) {
      ans += (SolvePhi(a / x, i + 1) + prime[i]) * f;
    }
  }
  return ans;
}
int SolveMu(unsigned a, int b) {
  if(a < prime[b]) {
    return 0;
  }
  int ans = h[Id(a)] + (b - 1);
  for(unsigned i = b; i <= cnt && prime[i]*prime[i] <= a; ++i) {
    ans -= SolveMu(a / prime[i], i + 1);
  }
  return ans;
}
int main() {
  scanf("%d", &T);
  while(T--) {
    scanf("%u", &n), sn = sqrt(n);
    if(!n) {
      puts("0 0");
      continue;
    }
    cnt = id = 0;
    for(unsigned i = 1; i <= n; i = a[id] + 1) {
      a[++id] = n / (n / i), g[id] = (ll)a[id] * (a[id] + 1) / 2 - 1, h[id] = a[id] - 1;
    }
    for(unsigned i = 2; i <= sn; ++i) {
      if(h[i] != h[i - 1]) {
        prime[++cnt] = i, sum[cnt] = sum[cnt - 1] + i;
        unsigned sq = i * i;
        for(int j = id; a[j] >= sq; --j) {
          unsigned t = Id(a[j] / i);
          g[j] -= i * (g[t] - sum[cnt - 1]);
          h[j] -= h[t] - (cnt - 1);
        }
      }
    }
    for(int i = 1; i <= id; ++i) {
      g[i] -= h[i], h[i] *= -1;
    }
    // 上面的计算都是不考虑 1 的函数值的，要加上去
    printf("%lld %d\n", SolvePhi(n, 1) + 1, SolveMu(n, 1) + 1);
  }
  return 0;
}
```

### 参考来源
[Min-25筛学习笔记](https://cmxrynp.github.io/2018/12/03/Min-25%E7%AD%9B%E5%AD%A6%E4%B9%A0%E7%AC%94%E8%AE%B0/)