# String

## Suffix Array

* usage:

```cpp
    // s原数组, sa, n数组长度, m值域[0, m)
    SuffixArray::da(s, sa, n, m);
    // SuffixArray::rank, SuffixArray::height
    calheight(r, sa, n);
```

* sa数组: $sa_i$为$Suffix_i$在字符串中第一次出现的位置

* height数组: $height_i = LongestPrefix(suffix(sa_{i-1}), suffix(sa_{i}))$

* $suffix_j$和$suffix_k$的最长公共前缀为$min(height_{rank_j+1}, height_{rank_j+2}, ..., height_{rank_k})$

* 求字符串中某两个后缀的最长公共前缀: 可转化为求height中某区间的最小值， 即RMQ问题

* 可重叠最长重复子串: $ans = \max {height_i}$

* 不可重叠最长重复子串: 二分答案，扫height数组，将连续的满足$height_i \geq k$的后缀分为一组，若同一组中存在$sa_{max} - sa_{min} \geq k$, 说明存在一组不重叠的重复子串

* 可重叠的k次最长重复子串: 二分答案，扫height数组，分组，判断是否存在个数不小于k的组

* (待续)
