/*
    Adapting Simpson
    求形如：
    /a
    |  F(x)dx
    /b
    的积分
*/

namespace AdaptingSimpson {
  template<typename Functor>
  double simpson(const double &a, const double &b, const Functor &functor) {
    double c = (a + b) / 2;
    return (functor(a) + 4 * functor(c) + functor(b)) * (b - a) / 6;
  }
  template<typename Functor>
  double asr(double a, double b, double eps, double A, const Functor &functor) {
    double c = (a + b) / 2;
    double L = simpson(a, c, functor), R = simpson(c, b, functor);
    if (fabs(L + R - A) <= 15 * eps) {
      return L + R + (L + R - A) / 15.0;
    }
    return asr(a, c, eps / 2, L, functor) + asr(c, b, eps / 2, R, functor);
  }
  template<typename Functor>
  double integrate(const Functor &functor, double from, double to, double eps=1e-5) {
    return asr(from, to, eps, simpson(from, to, functor), functor);
  }
}
using AdaptingSimpson::integrate;

/*
使用举例：
double f(double x) {
  return sin(x);
}
int main() {
  cout << fixed << setprecision(1) << integrate(f, 0, 1) << '\n';
}
*/