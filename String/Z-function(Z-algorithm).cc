// source: https://cp-algorithms.com/string/z-function.html

#include <algorithm>
#include <vector>
#include <string>
using namespace std;

// 与exkmp类似
// z[i]: suf_s[i..n-1] 与 s[0...n-1] 的公共前缀长度
vector<int> getZ(string s) {
  int n = (int)s.length();
  vector<int> z(n);
  for (int i = 1, l = 0, r = 0; i < n; i++) {
    if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) ++z[i];
    if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
  }
  return z;
}