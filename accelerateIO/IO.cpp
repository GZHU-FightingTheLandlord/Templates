#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;

namespace SpeedUp {
#define BUFLEN 100000
	inline char ch() {
		static char buf[BUFLEN], *st = nullptr, *ed = nullptr;
		if (st == ed) {
			st = buf;
			ed = st + fread(buf, 1, BUFLEN, stdin);
			if (st == ed) {
				return -1;
			}
		}
		return *st++;
	}
	inline bool blank(char& x) {
		return (x == '\n' || x == ' ' || x == '\r' || x == '\t');
	}
   	inline bool read(int& x) {
		char c;
		while (blank(c = ch()));
		if (c == -1) return false;
		for (x = c - '0'; (c = ch()) >= '0' && c <= '9'; x = x * 10 + c - '0');
		return true;
	}
#undef BUFLEN
}
using SpeedUp::read;
