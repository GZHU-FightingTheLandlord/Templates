namespace rpn {
  typedef long long rpn_type;
  struct node {
    int what;
    rpn_type val;
    node(int _what = 0, rpn_type _val = 0) : what(_what), val(_val) {}
  };
  rpn_type calc(int what, const rpn_type &a, const rpn_type &b) {
    if(what == 1) return a + b;
    if(what == 2) return a - b;
    if(what == 3) return a * b;
    if(what == 4) return a / b;
    assert(false);
  }
  const int LeftBucket = -114514, RightBucket = -1919810;
  int toWhat(const char c) {
    switch(c) {
      case '(': return LeftBucket;
      case ')': return RightBucket;
      case '+': return 1;
      case '-': return 2;
      case '*': return 3;
      case '/': return 4;
    }
    assert(false);
  }
  rpn_type solve(const string &str) {
    const int n = str.size();
    queue<node> Queue;
    stack<int> Stack;
    for(int i = 0; i < n; i++) {
      if(isdigit(str[i])) {
        rpn_type cur = 0;
        int j = i;
        while(isdigit(str[j])) {
          cur = cur * 10 + str[j] - '0';
          j++;
        }
        i = j - 1;
        Queue.push(node(0, cur));
      } else {
        const int cur = toWhat(str[i]);
        if(cur == LeftBucket) {
          Stack.push(cur);
        } else if(cur == RightBucket) {
          while(Stack.top() != LeftBucket) {
            Queue.push(node(Stack.top()));
            Stack.pop();
          }
          Stack.pop();
        } else {
          while(!Stack.empty() && Stack.top() >= cur) {
            Queue.push(node(Stack.top()));
            Stack.pop();
          }
          Stack.push(cur);
        }
      }
    }
    while(!Stack.empty()) {
      Queue.push(node(Stack.top()));
      Stack.pop();
    }
    stack<rpn_type> forCalc;
    while(!Queue.empty()) {
      const node cur = Queue.front();
      Queue.pop();
      if(cur.what) {
        const rpn_type b = forCalc.top(); forCalc.pop();
        const rpn_type a = forCalc.top(); forCalc.pop();
        forCalc.push(calc(cur.what, a, b));
      } else {
        forCalc.push(cur.val);
      }
    }
    return forCalc.top();
  }
}
