/*
    Author: SemonChan
    Time  : 2018 / 4 / 16 00:06
    Remark：由于减法未开发， 此模板暂时只对无符号大整数有效。
*/

#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;

struct BigInt{
    static const int MAX = 5000;
    int len, neg;
    int v[MAX];
    BigInt() { len = neg = 0; memset(v, 0, sizeof v);   }
};

int read(BigInt& a);
BigInt plu(const BigInt& a, const BigInt& b);
BigInt red(const BigInt& a, const BigInt& b);
void println(const BigInt& v);

int read(BigInt& a)
{
    a = BigInt();
    int len, i = 0, j;
    char v[BigInt::MAX];
    if (scanf("%s", v) == EOF)
        return EOF;
    if (v[0] == '-')
        a.neg = 1, i = 1;
    a.len = len = strlen(v) - a.neg;
    for (j = len - 1; j >= 0; i++, j--)
        a.v[j] = v[i] - '0';
    return 1;
}

BigInt plu(const BigInt& a, const BigInt& b)
{
    BigInt c;
    if ((a.neg && b.neg) || (!a.neg && !b.neg))
    {
        c.neg = a.neg;
        c.len = max(a.len, b.len);
        for (int i = 0; i < c.len; i++)
        {
            c.v[i] += a.v[i] + b.v[i];
            c.v[i + 1] = c.v[i] / 10;
            c.v[i] %= 10;
        }
        while (c.v[c.len] > 0)
        {
            c.v[c.len + 1] = c.v[c.len] / 10;
            c.v[c.len] %= 10;
            c.len++;
        }
    }
    else if (a.neg && !b.neg)
    {
        c = a;
        c.neg = 0;
        return red(b, c);
    }
    else
    {
        c = b;
        c.neg = 0;
        return red(a, c);
    }
    return c;
}

/*  减法待开发       */
BigInt red(const BigInt& a, const BigInt& b)
{
    return a;
}

void println(const BigInt& v)
{
    if (v.neg)
        printf("-");
    for (int i = v.len - 1; i >= 0; i--)
        printf("%d", v.v[i]);
    printf("\n");
}

int main()
{
    BigInt a, b;
    while (read(a) != EOF)
    {
        read(b);
        println(plu(a, b));
    }
    return 0;
}

/*
    clang++ -std=c++11 calc.cpp
    ./a.out
*/
