#include<bits/extc++.h>
// 下面是set的例子 map的话将__gnu_pbds::null_type改成想要的结构即可
typedef __gnu_pbds::tree<int, __gnu_pbds::null_type, less<int>, __gnu_pbds::rb_tree_tag, __gnu_pbds::tree_order_statistics_node_update> ordered_set;
/* 
  iterator find_by_order(size_type order)
    找到第order+1小的迭代器，如果order太大会返回end()
  size_type order_of_key(const_key_reference r_key)
    询问这个tree中有多少个比r_key小的元素
  void join(tree &other)
    把other中所有元素移动到*this中（要求原来other和*this的key不能相交，否则会抛出异常）
*/
typedef __gnu_pbds::priority_queue<int, less<int>> pq;
/*
  void join(priority_queue &other)
    把other合并到*this，然后other会被清空
*/
