/*
链接：http://acm.nyist.edu.cn/JudgeOnline/problem.php?pid=35
*/

#include<bits/stdc++.h>
using namespace std;

const int MXLEN = 1000 + 5;
int fst[MXLEN];
char str[MXLEN];

typedef double myClass;
typedef myClass CSS;
// 需要给myClass重载各种运算符

CSS jud(int begin, int end)
{
    /*规定区间[begin, end]的优先级标准为fst[begin]*/
    int i;
    CSS k;
    for (i = begin; i <= end; i++)
    {
        if (str[i] == '+' && fst[i] == fst[begin])
        {
            k = jud(begin, i - 1) + jud(i + 1, end);
            return k;
        }
    }
    for (i = end; i >= begin; i--)
    {
        if (str[i] == '-' && fst[i] == fst[begin])
        {
            k = jud(begin, i - 1) - jud(i + 1, end);
            return k;
        }
    }
    for (i = begin; i <= end; i++)
    {
        if (str[i] == '*' && fst[i] == fst[begin])
        {
            k = jud(begin, i - 1) * jud(i + 1, end);
            return k;
        }
    }
    for (i = end; i >= begin; i--)
    {
        if (str[i] == '/' && fst[i] == fst[begin])
        {
            k = jud(begin, i - 1) / jud(i + 1, end);
            return k;
        }
    }
    /*
    按照上面优先级： '/' > '*' > '-' > '+'
    添加第5种运算：
        for (i = end; i >= begin; i--)
        {
            if (str[i] == '^' && fst[i] == fst[begin])
            {
                k = jud(begin, i - 1) ^ jud(i + 1, end);
                return k;
            }
        }
    */
    if (str[begin] == '(')
    {
        for (i = begin + 1; fst[i] >= fst[begin + 1]; i++);

        k = jud(begin + 1, i - 1);
    }
    else
    {
        char *p = str;
        sscanf(p + begin, "%lf", &k);
    }
    return k;
}

CSS solve()
{
    const int len = strlen(str);
    for (int i = 1; i <= len - 1; i++)
    {
        if (str[i - 1] == '(')
        {
            fst[i] = fst[i - 1] + 1;
        }
        else
        {
            if (str[i] == ')')
            {
                fst[i] = fst[i - 1] - 1;
            }
            else
            {
                fst[i] = fst[i - 1];
            }
        }
    }
    return jud(0, len);
}

int main()
{
    int T;
    scanf("%d", &T);
    while (T--)
    {
        scanf("%s", str);
        const int len = strlen(str);
        str[len - 1] = '\0';
        myClass x = solve();
        printf("%.2lf\n", x);
    }
    return 0;
}
