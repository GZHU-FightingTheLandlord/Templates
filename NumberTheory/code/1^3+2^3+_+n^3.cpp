/*
  根据数学归纳法有:
  1^3 + 2^3 + 3^3 +...+ n^3 = (n * (n + 1) / 2) ^ 2
*/
#include <stdio.h>
int main()
{
  int n; scanf("%d", &n);
  int ans = n * (n + 1) / 2;
  printf("%d\n", ans * ans);
  return 0;
}
