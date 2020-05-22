// 最小表示法
int minPresentation(const string &s) {
  const int n = s.size();
  int i = 0, j = 1, k = 0;
  while (k < n && i < n && j < n) {
    if (s[(i + k) % n] == s[(j + k) % n]) {
      ++k;
    } else {
      s[(i + k) % n] > s[(j + k) % n] ? (i += k + 1) : (j += k + 1);
      i += (i == j);
      k = 0;
    }
  }
  return min(i, j);
}
