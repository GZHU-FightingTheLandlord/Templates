pll exgcd(const long long x, const long long y) {
  if (!y) return {1, 0};
  pll cur = exgcd(y, x % y);
  return {cur.second, cur.first - (x / y) * cur.second};
}
