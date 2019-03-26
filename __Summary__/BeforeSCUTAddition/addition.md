## 莫比乌斯函数线性筛

```cpp
    const int MAXN = 1e7;
    bool isnpri[MAXN];
    vector<int>pri;
    int mu[MAXN];
    void mobius() {
        mu[1] = 1;
        for(int i = 2; i < MAXN; i++) {
            if(!isnpri[i]) {
                pri.push_back(i);
                mu[i] = -1;
            }
            for(int j = 0; j < int(pri.size()); j++) {
                const int cur = i * pri[j];
                if(cur >= MAXN) {
                    break;
                }
                isnpri[cur] = true;
                if(i % pri[j] == 0) {
                    mu[cur] = 0;
                    break;
                } else {
                    mu[cur] = -mu[i];
                }
            }
        }
    }
```

