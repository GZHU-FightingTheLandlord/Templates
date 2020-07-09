int tot, fail[N], son[N][26];

int encode(char c) {
  return int(c - 'a');
}
void initNode(int at) {
  fail[at] = 0;
  memset(son[at], 0, sizeof son[at]);
}
void init() {
  tot = 0, initNode(0);
}
int insert(const string &s) {
  int cur = 0;
  for (auto i : s) {
    int c = encode(i);
    if (!son[cur][c]) {
      son[cur][c] = ++tot;
      initNode(tot);
    }
    cur = son[cur][c];
  }
  return cur;
}
void build() {
  queue<int> Q;
  for (int i = 0; i < 26; i++) {
    if (son[0][i]) {
      Q.push(son[0][i]);
    }
  }
  while (!Q.empty()) {
    int u = Q.front();
    Q.pop();
    for (int i = 0; i < 26; i++) {
      if (son[u][i]) {
        fail[son[u][i]] = son[fail[u]][i];
        Q.push(son[u][i]);
      } else {
        son[u][i] = son[fail[u]][i];
      }
    }
  }
}