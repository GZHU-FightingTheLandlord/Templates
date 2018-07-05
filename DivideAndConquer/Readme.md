# 分治法

将某一较大问题分割为一系列较小问题， 并通过解决较小问题得到原问题的解

## 二分搜索

### 基本思想

在有序数组中查找某一元素时，每次取中间元素进行比较，根据比较结果判断下一步应在哪半边进行
查找。

### 常用模板

* std::lower_bound

```cpp
    // 返回指向大于等于value的第一个元素的迭代器
    template <class ForwardIt, class T>
    ForwardIt lower_bound(ForwardIt first, ForwardIt last, const T& value);
```

```cpp
    // comp返回值: true表示a应在b前面
    // comp函数接口:
    // bool comp(const TYPE1& a, const TYPE2& b);
    template <class ForwardIt, class T, class Compare>
    ForwardIt lower_bound(ForwardIt first, ForwardIt last, const T& value, Compare comp);
```

* std::upper_bound

```cpp
    // 返回指向大于value的第一个元素的迭代器
    template <class ForwardIt, class T>
    ForwardIt upper_bound(ForwardIt first, ForwardIt last, const T& value);
```

```cpp
    // comp函数参照std::lower_bound
    template <class ForwardIt, class T>
    ForwardIt upper_bound(ForwardIt first, ForwardIt last, const T& value, Compare comp);
```

## 最近点对

懒得写，直接贴链接:

* [Wiki](https://en.wikipedia.org/wiki/Closest_pair_of_points_problem)

* [GeeksForGeeks](https://www.geeksforgeeks.org/closest-pair-of-points/)
