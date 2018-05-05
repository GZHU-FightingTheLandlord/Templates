/*
  设 m = n ^ n
  对两边取log10有 log(m) = n * log(n)
  则 m = 10 ^ (n * log(n))
  n * log(n) 可记为 N + S, 其中N为整数部分, S为小数部分
  因10的整数次方首位必为1, 故首位数字与小数部分有关
  即 ans = 10 ^ S
*/
#include <stdio.h>
#include <math.h>
int main()
{
  double n; scanf("%lf", &n);
  double e = n * log10(n);
  e = fmod(e , 1.0);
  printf("%d\n", int(pow(10, e));
  return 0;
}
