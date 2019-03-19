/*
    Adapting Simpson
    求形如：
    /a
    |  F(x)dx
    /b
    的积分
*/

/* 
   调用： asr(a, b, eps, simpson(a, b))
   其中eps为自定义误差
*/ 

double simpson(const double& a, const double& b)
{
    double c = (a + b) / 2;
    return (F(a) + 4 * F(c) + F(b)) * (b - a) / 6;
}

double asr(double a, double b, double eps, double A)
{
    double c = (a + b) / 2;
    double L = simpson(a, c), R = simpson(c, b);
    if (fabs(L + R - A) <= 15 * eps)
        return L + R + (L + R - A) / 15.0;
    return asr(a, c, eps / 2, L) + asr(c, b, eps / 2, R);
}
