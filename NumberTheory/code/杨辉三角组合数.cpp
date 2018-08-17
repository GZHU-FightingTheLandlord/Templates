#define Rep(i, a, b) for(int i = (a); i <= (b); i++)

const int MOD = 1e9 + 7;

int C[1005][1005];

void init() {
	Rep(i, C[0][0] = 1, 1004) Rep(j, C[i][0] = 1, i) {
		C[i][j] = (C[i - 1][j - 1] + C[i - 1][j]) % MOD;
	}
}