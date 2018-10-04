#include <string.h>

const int maxn = 1e6 + 5;

int fail[maxn];

void getfail(char *par, int n) {
	for (int i = 0, j = fail[0] = -1; i < n; i++, j++) {
		while (~j && par[j] != par[i]) j = fail[j];
		fail[i + 1] = j + 1;
	}
}

int solve(char *str, int n, char *par, int m) {
	for (int i = 0, j = 0; i < n; ) {
		while (~j && str[i] != par[j]) j = fail[j];
		++i, ++j;
		if (j >= m) return i - m + 1;
	}
	return -1;
}