# NumberTheory

## 九余数定理

* 一个数各位数字之和等于这个数对9取模所得的数

* 每次将指数进行一次$log(N)$级别的变换
* 矩阵快速幂: 在$O(log(N))$级别的时间求第n个斐波納契数列$f(n)=a*f(n-1)+b*f(n-2)$
* 快速乘: 利用二进制实现$ab \mod p$, 防止溢出

## 母函数(组合数学)

* 详见模板和实例

## 五边形数定理

* 五边形数: 1, 5, 12, 22, ……
* 第(n-1)个三角数+n^2为第n个五边形数
* $S[n] = s[n-1]+3n-2$

## zeckendorf定理

* 任何正整数可以表示为若干个不连续的Fibonacci数之和(斐波納契博弈)

## 错排公式

* $d(n) = (n - 1)(d[n-2] + d[n-1])$

## 不定方程

* 二元一次不定方程ax+by=c有解的充要条件是(a,b)|c

## 欧拉定理

* 欧拉函数: φ(n)是小于等于n的正整数中与n互质的数的数目
* 若n, a为正整数且n与a互质, 则a^φ(n)≡1(mod n)
* 费马小定理(Fermat小定理): 对任意a和任意质数p: a^p≡a(mod p), 若a不能被p整除, a^(p-1)≡1(mod p)
* 欧拉降幂: $$x^y\mod p = x^{y \mod \phi(p) \space + \space \phi(p)}, y \geq \phi(p)$$

```py
    # (x ^ y) % p
    def calc(x, y, p):
        if y < p:
            return qpow(x, y, p)
        else:
            # (x ^ y) % p = (x ^ (y % φ(p) + φ(p))) % p
            return qpow(x, y % phi[p] + phi[p], p)
```

## 费马大定理 && 费马方程

* $x^n + y^n = z^n$, 由费马大定理可知: 若$n≥2$且$n$为整数, 则不存在整数解$(x,y,z)(xyz≠0)$

## 同余方程（符号待更新)

## SG(Sprague-Grundy)函数

* 对于任意状态x, 它的SG函数值g(x)=mex{g(y)|y是x的后续状态}, mex是一个对于非负整数集合S的运算, mex(S)为S中没有出现的最小非负整数
* 终止状态的SG函数值为0
* 博弈打表(待更新)

## Lucas定理

* $C(n, m) \mod p = C(n / p, m / p) C(n \mod p, m \mod p) \mod p$, 其中p为质数

## pick定理(实际上属于计算几何)

* 给定正方形格子点的简单多边形, i为其内部格点数目, b为其边上格点数, 则其面积$A=i+b/2-1$

## 逆元(费马小定理)

* 当p为素数时， $a / b \mod p = ab_1 \mod p,  b_1 = b ^ {p - 2} \mod p$

## 威尔逊定理

* 当且仅当p为素数时: $(p-1)!≡-1(\mod p)$

## 杨辉三角(应用于排列组合)

* 第n行的元素个数有n个
* 第n行的元素之和为2^(n-1)
* 第n行第m个数的值为C(n-1,m-1),C为组合数
* (a+b)^n展开后的各项系数等于第n+1行的值
* 第n行第m个数的奇偶判断(m-1)&(n-1)==(m-1)?odd:even

## 第一类斯特林数

* $S(p, k)$表示把p个人分成k组作环排列的方案数

* $ S(p, k) = (p - 1) * S(p - 1, k) + S(p - 1, k - 1), 1 \leq k \leq p - 1$

* $S(p,0)=0, p\geq1$

* $S(p,p)=1, p\geq0$

* 使第p个人单独构成一个环排列，前p-1人构成k-1个环排列，方案数$S(p-1,k-1)$。

* 使前p-1个人构成k个环排列，第p个人插入到第i人左边，方案数$S(p-1, k)$

## 第二类斯特林数

* 将p个物体分成k个非空的不可辨别的集合的方案数

* $S(p, k) = k * S(p-1, k) + S(p-1, k-1), 1 \leq k \leq p - 1$

* $S(p,0)=0, p\geq1$

* $S(p,p)=1, p\geq0$

* 考虑第p个物品，p可以单独构成一个非空集合，此时前p-1个物品构成k-1个非空的不可辨别的集合，方法数为$S(p-1,k-1)$

* 也可以前p-1种物品构成k个非空的不可辨别的集合，第p个物品放入任意一个中，这样有$kS(p-1,k)$种方法。

## 拉格朗日插值公式

* $L(x):=\sum _{{j=0}}^{{k}}y_{j}\ell _{j}(x)$

* $\ell _{j}(x):=\prod _{{i=0,\,i\neq j}}^{{k}}{\frac  {x-x_{i}}{x_{j}-x_{i}}}={\frac  {(x-x_{0})}{(x_{j}-x_{0})}}\cdots {\frac  {(x-x_{{j-1}})}{(x_{j}-x_{{j-1}})}}{\frac  {(x-x_{{j+1}})}{(x_{j}-x_{{j+1}})}}\cdots {\frac  {(x-x_{{k}})}{(x_{j}-x_{{k}})}}.$

## 莫比乌斯反演公式

* $F(n)=\sum\limits_{d\mid n}f(d)\Rightarrow f(n)=\sum\limits_{d\mid n}\mu(d)F(\frac{n}{d})$

* $F(n)=\sum\limits_{n\mid d}f(d)\Rightarrow f(n)=\sum\limits_{n\mid d}\mu(\frac{d}{n})F(d)$
