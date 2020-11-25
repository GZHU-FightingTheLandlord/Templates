std::vector<int> getFail(const std::string &s) {
  const int n = (int) s.length();
  std::vector<int> fail(n + 1);
  for (int i = 0, j = fail[0] = -1; i < n; i++) {
    while (j != -1 and s[j] != s[i]) {
      j = fail[j];
    }
    fail[i + 1] = j + 1;
  }
  return fail;
}

int match(const std::string& s, const std::string& t, const std::vector<int>& fail) {
  const int n = (int) s.length(), m = (int) t.length();
  for (int i = 0, j = 0; i < n; i++, j++) {
    while (j != -1 and s[i] != t[j]) {
      j = fail[j];
    }
    if (j + 1 == m) {
      return (i + 1) + 1 - m;
    }
  }
  return -1;
}