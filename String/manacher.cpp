#include <string>
#include <vector>
#include <utility>

int Manacher(const std::string& t) {
  std::string s = "$#";
  for (auto c : t) {
    s += c, s += "#";
  }
  const int n = (int) s.length();
  int ans = 0, ind = 0, right = 0;
  std::vector<int> dp(n, 1);
  for (int i = 1; i < n; i++) {
    if (i < right) {
      dp[i] = std::min(right - i, dp[2 * ind - i]);
    }
    while (i + dp[i] < n && s[i + dp[i]] == s[i - dp[i]]) {
      ++dp[i];
    }
    ans = std::max(ans, dp[i] - 1);
    if (i + dp[i] > right) {
      ind = i, right = i + dp[i];
    }
  }
  return ans;
}