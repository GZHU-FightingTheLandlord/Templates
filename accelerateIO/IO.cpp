#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;

namespace SpeedUp {
#define BUFLEN 100000
	char buf[BUFLEN]; // Buffer
	size_t pos, end;

    // Get char from buffer
	inline char next()
	{
		if (pos == end) {
			pos = 0;
			end = fread(buf, 1, BUFLEN, stdin);
			if (pos == end) throw int(0); // throw EOFError (try catch)
		}
		return buf[pos++];
	}

    // Positive Interger
    template <typename T> inline void read(T& x)
	{
		char c;
		while ((c = next()) == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '-');
		x = c - 48;
		while ((c = next()) > 47 && c < 58) {
			x = x * 10 + c - 48;
		}
	}
#undef BUFLEN
}
using SpeedUp::read;