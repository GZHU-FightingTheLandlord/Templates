# BalaBala

## nth_element

    #include <algorithm>

    using std::nth_element;

### 基本模式/重载小于号

* 函数原型:

        template \<class Iter>

        void nth_eleent(Iter first, Iter nth,
        Iter last);

* 将序列中第n大元素移至nth指向位置， 且first到nth间元素均不大于\*nth， nth到last间元素均不小于\*nth

* 调用样例:

    int n = 5, k = 3;

    int v[5];

    nth_element(v, v + k - 1, v + n);

    printf("%d\n", v[2]); // 此时v数组中第3大元素移至v[2]中

### cmp函数重载

* 函数原型:

    template /<class Iter, class Compare>

    void nth_element(Iter first, Iter nth, Iter last, Compare Comp);

* cmp函数定义:

    // 返回true表示a应在b前面, false反之

    template \<typename T>

    bool cmp(const T& a, const T& b);

## MergeSort

### 工作原理<del>（瞎BB）</del>

递归地将原序列划分为多个小序列，然后将有序的小序列合并即可

### Example

* 原序列为: [5, 3, 1, 2, 4]

* 划分: <br /> [5, 3, 1, 2, 4] <br /> [5, 3], [1, 2], [4] <br /> [5], [3], [1], [2], [4]

* 合并: <br /> [3, 5], [1, 2], [4] <br /> [1, 2, 3, 5], [4] <br /> [1, 2, 3, 4, 5]

## QuickSort

### <del> 依然是瞎BB </del>

以从序列中选取的元素为标杆，将序列分为两部分（大于标杆或小于标杆）， 然后对划分出的两个序列重复这一操作， 直至有序

### Example

* 原序列: [5, 3, 1, 2, 4]

* 以第一个元素为标杆示例: <br /> [3, 1, 2, 4], 5 <br /> [1, 2], 3, [4], 5 <br /> 1, 2, 3, 4, 5
