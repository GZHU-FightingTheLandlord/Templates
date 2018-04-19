class Kmp
{
	int pre[MAX];
	char s[MAX], p[MAN];

	void init()
	{
		fill(pre, pre + MAX, -1);
		for(int i = 1, j = -1; p[i]; ++i)
		{
			while(j >= 0 && p[i] != p[j + 1])
				j = pre[j];
			if(p[i] == p[j + 1])
				++j;
			pre[i] = j;
		}
	}

	void solve(vector<int> &match)
	{
		match.clear(); prepare();
		for (int i = 0, j = -1; s[i]; ++i)
		{
			while(j >= 0 && s[i] != p[j + 1])
				j = pre[j];
			if(s[i] == p[j + 1])
				++j;
			if(!p[j + 1])
				match.push_back(i - j);
		}
	}
};
