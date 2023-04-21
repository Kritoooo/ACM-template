# 数论

## 最大公约数

### gcd

```c
int gcd(int a, int b) { return b ? gcd(b, a % b) : a; }
```

### lcm

```c
long long lcm(int a, int b) { return 1ll * a / gcd(a, b) * b; }
```

## 线性筛

​	从小到大枚举因子

​	p[i] : i的最小素因子

​	prime[i]：素数的值

```c
void init(int n) {
    p[1] = 1;
    for(int i = 2; i <= n; i ++) {
        if(!p[i]) p[i] = i, prime[++tot] = i;
        for (int j = 1; j <= tot && prime[j] * i <= n; j ++) {
            p[i * prime[j]] = prime[j];
            if(p[i] == prime[j]) break;
        }
    }
}
```



## 扩展欧几里得

​	用来求解不定方程
$$
ax+by=gcd(a,b)
$$

```c
int exgcd(int a, int b, int &x, int &y) {
    if(b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int d = exgcd(b, a % b, y, x);
    y -= (a / b) * x;
    return d;
}
```

​	求解非负整数解(x, y)，输出x最小的解
$$
ax+by=d
$$

```c
cin >> a >> b >> m;
ll d = exgcd(a, b, x, y);
if(m % d != 0) {
    cout << -1 << "\n";
    continue;
}
a /= d; b /= d; m /= d;
__int128 xx = (__int128)x * m;
xx %= b;
if (xx < 0) xx += b;
__int128 yy = (m - a * xx) / b;
if(yy < 0) {
    cout << -1 << "\n";
    continue;
}
cout << (ll)xx << " " << (ll)yy << "\n";
```



## 算数基本定理

​	任何一个大于1的自然数N，如果N不为质数，那么N可以唯一分解成有限个质数的乘积

​	
$$
N={p_{1}}^{e_{1}}*{p_{2}}^{e_{2}}*{p_{3}}^{e_{3}}...*{p_{n}}^{e_{n}} (p_{1} < p_{2} < p_{3} < ...<p_{n})
$$


## 欧拉函数

​	对正整数n，欧拉函数是小于n的正整数中与n互质的数的数目。
$$
\varphi (x) = x*\prod_{i = 1}^{n}(1-\frac{1}{p_{i}})
$$


```c
int phi(int n) {
    int phin = n, res = n;
    for(int i = 2;  i * i <= n; i++) {
        if(n % i == 0) {
            phin = phin / i * (i - 1);
            while(n % i == 0) n /= i;
        }
    }
    if(n > 1) phin = phin / n * (n - 1);
    return phin;
}
```



## 欧拉定理

设a,m都为正整数，且gcd(a,m) = 1,则有
$$
a^{\varphi (m)}\equiv 1(mod\ m)
$$


## 逆元

​	任意整数a和其逆元满足
$$
aa^{-1}\equiv 1(mod\ n)
$$
​	在对除法运算进行取模时，利用逆元有 (a / b) mod p = a * inv(b) (mod p).

​	当a为b的倍数时，有
$$
a÷bmodp=(amod(b×p))÷b
$$


​	费马小定理求逆元（模数为素数时）

```c
long long quickpow(long long a, long long b) {
    if (b < 0) return 0;
    long long ret = 1;
    a %= mod;
    while(b) {
        if (b & 1) ret = (ret * a) % mod;
        b >>= 1;
        a = (a * a) % mod;
    }
    return ret;
}
long long inv(long long a) {
    return quickpow(a, mod - 2);
}
```

​	扩展欧几里得求逆元

```c
int getInv(int a, int mod) {
    int x, y;
    int d = exgcd(a, mod, x, y);
    return d == 1 ? (x % mod + mod) % mod : -1;
}
```

​	欧拉定理同样可以求逆元

​	求解1到n的逆元

```c
    inv[1] = 1;
    for(int i = 2; i <= n; i++) {
        inv[i] = (p - p / i) * inv[p % i] % p;
    }
```

​	求解阶乘逆元

```c
void init()
{
    fac[0] = 1;
    for (int i = 1; i < maxn; i++)
    {
        fac[i] = fac[i - 1] * i % mod;
    }
    inv[maxn - 1] = quick_pow(fac[maxn - 1],mod - 2,mod);
    for (int i = maxn - 2; i >= 0; --i)
    {
        inv[i] = inv[i + 1] * (i + 1) % mod;
    }
}
```

求n个数的逆元

```c
    s[0] = 1;
    for (int i = 1; i <= n; i++) s[i] = s[i - 1] * a[i] % p;
    int x, y;
    exgcd(s[n], p, x, y);
    if (x < 0) x += p;
    t[n] = x;
    assert(s[n] * x % p == 1);
    for(int i = n; i >= 1; i --) t[i - 1] = t[i] * a[i] %p;
    for(int i = 1; i <= n; i++) {
        inv[i] = s[i - 1] * t[i] % p;
    }
```



## 中国剩余定理

## 整除分块

```c
for(ll l = 1; l <= n; l ++) {
    ll d = n / l, r = n / d;
    sum += (r - l + 1) * d;
    // l .. r  n / x = d
    l = r;
}
```

## long long 取模

```c
ll mul(ll x, ll y, ll m) {
    x %= m; y %= m;
    ll d = ((long double) x * y / m);
    d = x * y - d * m;
    if (d >= m) d -= m;
    if(d < 0) d += m;
    return d;
}
```



## Lucas定理

```c
long long Lucas(long long n, long long m, long long p) {
    if(m == 0) return 1;
    return (c(n % p, m % p, p) * Lucas(n / p, m / p, p)) % p;
}
```

## 扩展欧拉定理

$$
a^{b}\ \% \ m = a^{b \% \ \varphi(m) + \varphi(m)}\ \%\ m
$$



## 积性函数

定义：(a,b) = 1, f(ab) = f(a)f(b)

### 常见的积性函数

1.
$$
id(x) = x
$$
2.
$$
1(x) = 1
$$
3.
$$
e(x) = \left\{\begin{matrix}
 1& x = 1 & \\ 
 0& x\neq 1 & 
\end{matrix}\right.
$$
4.欧拉函数
$$
\varphi (n) = n\sum_{p | n}(1-\frac{1}{p})
$$
5.d(n)因子个数
$$
d(p^{e}) = e + 1
$$
6.
$$
\sigma  (p^{e}) = p^{0} + p^{1} + ...+p^{e}
$$
7.莫比乌斯函数
$$
\mu (p^{e})=\left\{\begin{matrix}
1 & e=0 \\ 
-1 & e=1\\ 
 0&e\geq 2 
\end{matrix}\right.
$$

性质
$$
\sum_{d|n}\mu(d)=[n=1]
$$


### 线性筛求积性函数

由
$$
f(n) = f(p_{1}^{e_{1}})*f(p_{2}^{e_{2}}) *...*f(p_{k}^{e_{k}})
$$
得
$$
f(n) = f(p_{1}^{e_{1}})*f(n/p_{1}^{e_{1}})
$$

```c
const int N = 2e7 + 1000;
int p[N], pr[N / 5], n, pe[N], tot;
int f[N], a, b, ans;
void prime() {
    p[1] = 1;
    for(int i = 2; i <= n; i ++) {
        if(!p[i]) p[i] = i, pe[i] = i, pr[++tot] = i;
        for(int j = 1; j <= tot && pr[j] * i <= n; j ++) {
            p[i * pr[j]] = pr[j];
            if (p[i] == pr[j]) {
                pe[i * pr[j]] = pe[i] * pr[j];
                break;
            } else {
                pe[i * pr[j]] = pr[j];
            }
        }
    }
}

void compute(int n, function<void(int)> calcpe) {
    f[1] = 1;
    for(int i = 2; i <= n; i ++) {
        if(i == pe[i]) calcpe(i);
        else f[i] = f[pe[i]] * f[i / pe[i]];
    }
}
//因子个数
compute(n, [&](int x){
        f[x] = f[x / p[x]] + 1;
});
//因子和
compute(n, [&](int x){
       f[x] = f[x / p[x]] + x;
});
//欧拉函数
compute(n, [&](int x){
      f[x] = x / p[x] * (p[x] - 1);	
});
//mo'b
compute(n, [&](int x){
      f[x] = x == p[x] ? -1 : 0;
});
```

### 不同质因子个数

​	求一个数不同的质因子个数 p[i]。（加性函数）

```c
void init(int n) {
    p[1] = 1;
    for(int i = 2; i <= n; i ++) {
        if(!vis1[i]) p[i] = 1, prime[++tot] = i;
        for (int j = 1; j <= tot && prime[j] * i <= n; j ++) {
            vis1[i * prime[j]] = 1;
            if(i % prime[j] == 0) {
                p[i * prime[j]] = p[i];
                break;
            }
            p[i * prime[j]] = p[i] + p[prime[j]];
        }
    }
}
```



## 迪利克雷卷积

### 定义

$$
h(n) = \sum_{d|n}f(d)g(n/d) = \sum_{d_{1}d_{2} = n}f(d_{1})g(d_{2})
$$

常见的卷积

1.d(n) 因子个数
$$
d(n) = \sum_{d|n}1(d)1(d/n)
$$
2.因子和
$$
\sigma(n) = \sum_{d|n}id(d)1(n/d)
$$
3.
$$
f = f*e
$$
4.
$$
e = 1*\mu
$$


### 性质

1.交换律
$$
h = f*g = g*f
$$
2.结合律
$$
p=(f*g)*h=f*(g*h)
$$
**3.f和g是积性函数，则f*g也是积性函数**

## 莫比乌斯反演

### 形式

$$
f(n) = \sum_{d|n}g(d) <=>g(n) = \sum \mu(n/d)f(d)
$$

$$
f = g*1 <=>g=f*u
$$

### 一些经典的反演

1.

![image-20221101224706988](C:\Users\xsf\AppData\Roaming\Typora\typora-user-images\image-20221101224706988.png)

在求
$$
\sum^{n}_{i=1}\sum^{m}_{j=1}[gcd(i,j)=1](n<m)
$$
时，可以通过上式变化为![image-20221101225037605](C:\Users\xsf\AppData\Roaming\Typora\typora-user-images\image-20221101225037605.png)化简得![image-20221101225209016](C:\Users\xsf\AppData\Roaming\Typora\typora-user-images\image-20221101225209016.png)该式可以通过整除分块O(根号n)求解

形如![image-20221101225357216](C:\Users\xsf\AppData\Roaming\Typora\typora-user-images\image-20221101225357216.png)同样可以通过上述方法求解

code如下：

```c
//P3455
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const int N = 1e6 + 1000;
int p[N], pr[N / 5], n, pe[N], tot;
int f[N], smu[N], a, b, ans;
void prime() {
    p[1] = 1;
    for(int i = 2; i <= n; i ++) {
        if(!p[i]) p[i] = i, pe[i] = i, pr[++tot] = i;
        for(int j = 1; j <= tot && pr[j] * i <= n; j ++) {
            p[i * pr[j]] = pr[j];
            if (p[i] == pr[j]) {
                pe[i * pr[j]] = pe[i] * pr[j];
                break;
            } else {
                pe[i * pr[j]] = pr[j];
            }
        }
    }
}

void compute(int n, function<void(int)> calcpe) {
    f[1] = 1;
    for(int i = 2; i <= n; i ++) {
        if(i == pe[i]) calcpe(i);
        else f[i] = f[pe[i]] * f[i / pe[i]];
    }
}//莫比乌斯函数求解

int main() {
    int q;
    cin >> q;
    n = 5e4 + 10;
    prime();
    compute(n, [&](int x){
      f[x] = x == p[x] ? -1 : 0;
    }); 
    for(int i = 1; i <= n; i ++) {
        smu[i] = smu[i - 1] + f[i];
    }
    while(q --) {
        int a, b, k;
        cin >> a >> b >> k;
        a /= k, b /= k;
        if(a > b) swap(a, b);
        ll ans = 0;
        for(int l = 1; l <= a; l ++) {
            int r = min(a / (a / l), b / (b / l));
            ans += (1ll*smu[r] - smu[l - 1]) * (a / l) * (b / l);
            l = r; 
            // cout << l << endl;
        }//整除分块
        cout << ans << endl;
    }
}
```



## 高斯消元

### 异或高斯消元	

```c
int xor_gauss() {
    int r = 1;
    for (int c = n; c >= 1; c--) {
        int temp = r;
        for (int i = r; i <= n; i++)
            if (a[i][c]) {
                temp = i;
                break;
            }
        if (!a[temp][c]) continue;
        swap(a[r], a[temp]);
        for (int i = r + 1; i <= n; i++) {
            if (a[i][c]) a[i] ^= a[r];
        }
        r++;
    }
    if (r <= n) {
        for (int i = 1; i <= n; i++) {
            if (a[i] == 1) return -1;
            if (a[i] == 0) return quick_pow(2, n - i  + 1, mod);
        }
    }
    for (int i = 1; i <= n; i++) {
        int res = a[i][0];
        for (int j = i - 1; j; j--)
            res ^= a[i][j] ^ ans[j];
        ans[i] = res;
    }
    return 1;
}
```

### 高斯消元求行列式 （取模）

```c
int gauss(int a[N][N],int n)
{
	int res=1;
	for(int i=1; i<=n; i++)
	{
		for(int j=i+1; j<=n; j++)
		{
			while(a[j][i])
			{
				int t=a[i][i]/a[j][i];
				for(int k=i; k<=n; k++)
					a[i][k]=(a[i][k]-t*a[j][k]+mod)%mod;
				swap(a[i],a[j]);
				res=-res;
			}
		}
		res=(res*a[i][i])%mod;
	}
	return (res+mod)%mod;
}
```



# 图论

## 邻接表建图

```c
void add(int a, int b) {
    e[idx] = b, ne[idx] = h[a], h[a] = idx ++;
}
```

## 最短路

### floyd

#### 传递闭包

```c
for(int k = 1; k <= n; k ++)
        for(int i = 1; i <= n; i ++)
            for(int j = 1; j <= n; j ++)
                d[i][j] |= d[i][k] & d[k][j];
```



## LCA

### 倍增求LCA

```c
int fa[N][16]
void bfs() {
    memset(depth, 0x3f, sizeof depth);
    depth[root] = 1, depth[0] = 0;
    int hh = 0, tt = -1;
    q[++tt] = root;
    while(hh <= tt) {
        int t = q[hh ++];
        for(int i = h[t]; ~i; i = ne[i]) {
            int j = e[i];
            if(depth[j] > depth[t] + 1) {
                depth[j] = depth[t] + 1;
                q[++tt] = j;
                fa[j][0] = t;
                for(int k = 1; k <= 15; k++) 
                    fa[j][k] = fa[fa[j][k - 1]][k - 1];
            }
        }
    }
}

int lca(int a, int b) {
    if(depth[a] < depth[b]) swap(a, b);
    for(int k = 15; k >= 0; k--) 
        if(depth[fa[a][k]] >= depth[b]) 
            a = fa[a][k];
    if(a == b) return a;
    for(int k = 15; k >= 0; k --)
        if(fa[a][k] != fa[b][k]) {
            a = fa[a][k];
            b = fa[b][k];
        }
    return fa[a][0];
}
```

### 树上差分

点差分

将两点u,v之间路径上的所有点权增加x, o = LCA(u, v), o的父亲结点为p，则：
$$
diff[u]+=x,diff[v]+=x,diff[o]-=x,diff[p]-=x;
$$
code

```c
void dfs(int u, int fa) {  
    int res = 0;  
    for(int i = h[u]; ~i; i = ne[i]) {  
        int j = e[i];  
        if(j == fa) continue;  
        dfs(j, u);  
        res += sum[j];  
    }   
    sum[u] = res + diff[u];  
}    
cin >> u >> v;  
int x = lca(u, v);  
diff[u] += 1;  
diff[v] += 1;   
diff[x] -= 1;  
diff[p[x]] -= 1;  
```

边差分

将两点u, v之间路径上的边权增加x，o = LCA(	u,v)，以每条边两端深度较大的节点存储改边的差分数组，则操作如下：
$$
diff[u]+=x,diff[v]+=x,diff[o]-=2*x;
$$


## Tarjan缩点

### 有向图的强连通分量

连通分量：对于分量中任意两点u，v必然可以从u走到v，且从v走到u

强连通分量：极大连通分量

通过Tarjan缩点可以让有向图转变为有向无环图，转变后的图里面的每一个点是原图的一个强连通分量。

一个有向图，变成一个强连通分量至少需要添加max(p,q)条边，p为缩点后入度为0的点，q为缩点后出度为0的点

```c
int dfn[N], low[N], timestamp;
int stk[N], top;
bool in_stk[N];
int id[N], scc_cnt, size[N];

void tarjan(int u) {
    dfn[u] = low[u] = ++ timestamp;
    stk[ ++ top ] = u, in_stk[u] = 1;
    for(int i = h[u]; i != -1; i = ne[i]) {
        int j = e[i];
        if(!dfn[j]) {
            tarjan(j);
            low[u] = min(low[u], low[j]);
        }
        else if(in_stk[j]) low[u] = min(low[u], dfn[j]);
    }

    if(dfn[u] == low[u]) {
        ++ scc_cnt;
        int y;
        do {
            y = stk[top -- ];
            in_stk[y] = 0;
            id[y] = scc_cnt;
            size[scc_cnt] ++;
        } while(y != u);
    } 
}
```

### 差分约束

```c
for(int i = 1; i <= m; i++) {
    int f, x, y;
    std::cin >> f >> x >> y;
    if(f == 1) { // 等于 x == y, x -> y, w = 0, y -> x, w = 0;
        add(h,x,y,0);add(h,y,x,0);
    } else if(f == 2) { // 小于 x < y, x->y, w = 1;
        add(h,x,y,1);
    } else if(f == 3) { // 大于等于 x >= y, y->x, w = 0;
        add(h,y,x,0);
    } else if(f == 4) { // 大于 x > y, y->x, w = 1;
        add(h,y,x,1);
    } else if(f == 5) { // 小于等于 x <= y, x->y, w = 0;
        add(h,x,y,0);
    } 
}

for(int i = 1; i <= n; i++) {
    add(h,0,i,1);
}
// 缩点
tarjan(0);
// 建DAG
int f = 1;
for(int i = 0; i <= n; i ++) {
    for(int j = h[i]; ~j; j = ne[j]) {
        int k = e[j];
        int u = id[i], v = id[k];
        if(u == v) {
            // 不符合不等式关系
            if(w[j] > 0) {
                f = 0;
                break;
            }
        } else add(h1, u, v, w[j]);
    }
    if(!f) break;
}
if(!f) {
    std::cout << "-1\n";
    return 0;
}
// DAG求最长路
for(int i = scc_cnt; i >= 0; i --) {
    for(int j = h1[i]; ~j; j = ne[j]) {
        int k = e[j];
        d[k] = std::max(d[i] + w[j], d[k]);
    }
}
```

​	

### 无向图的双连通分量

#### 边双连通分量 e-Dcc

​	无向图中，极大的不含有桥的连通块被称为边的双连通分量

​	在里面不管删掉哪条边，仍然连通。

​	每对点之间至少存在两条没有公共边的路径

​	一个无向图变成边双连通分量，至少需要添加(cnt + 1)/ 2 条边，cnt为缩点后度数为1的点的个数

```c
int dfn[N], low[N], timestamp;  // 时间戳
int stk[N], top;
int id[N], dcc_cnt;  // 每个点所属分量编号
bool is_bridge[M];

void tarjan(int u, int from)
{
    dfn[u] = low[u] = ++ timestamp;
    stk[ ++ top] = u;
    
    for (int i = h[u]; ~i; i = ne[i])
    {
        int j = e[i];
        if (!dfn[j])
        {
            tarjan(j, i);
            low[u] = min(low[u], low[j]);
            if (dfn[u] < low[j])
                is_bridge[i] = is_bridge[i ^ 1] = true;
        }
        else if (i != (from ^ 1))
            low[u] = min(low[u], dfn[j]);
    }
    
    if (dfn[u] == low[u])
    {
        ++ dcc_cnt;
        int y;
        do {
            y = stk[top -- ];
            id[y] = dcc_cnt;
        } while (y != u);
    }
}
```



#### 点双连通分量 v-Dcc

​	极大的不包含割点的连通块被称为点的双连通分量

​	每个割点至少属于两个点双连通分量

```c
int h[N], e[M], ne[M], idx;
int dfn[N], low[N], timestamp;  // 时间戳
int stk[N], top;
int dcc_cnt;
vector<int> dcc[N];  // 每个分量有哪些点
bool cut[N];  // 是否为割点
int root;

void tarjan(int u)
{
    dfn[u] = low[u] = ++ timestamp;
    stk[ ++ top] = u;
    
    if (u == root && h[u] == -1)
    {
        dcc_cnt ++ ;
        dcc[dcc_cnt].push_back(u);
        return;
    }
    
    int cnt = 0;
    for (int i = h[u]; ~i; i = ne[i])
    {
        int j = e[i];
        if (!dfn[j])
        {
            tarjan(j);
            low[u] = min(low[u], low[j]);
            if (dfn[u] <= low[j])
            {
                cnt ++ ;
                if (u != root || cnt > 1) cut[u] = true;
                ++ dcc_cnt;
                int y;
                do {
                    y = stk[top -- ];
                    dcc[dcc_cnt].push_back(y);
                } while (y != j);
                dcc[dcc_cnt].push_back(u);
            }
        }
        else low[u] = min(low[u], dfn[j]);
    }
}

```



## 二分图

1.二分图不存在奇数环，染色法不存在矛盾

2.匈牙利算法，匹配，最大匹配，匹配点，增广路径

3.最大匹配数 = 最小点覆盖 = 总点数 - 最大独立集 = 总点数 - 最小路径覆盖

### 匈牙利算法

```c
bool find(int u) {  
    for(auto v : eg[u]) {  
        if(!st[v]) {  
            st[v] = 1;  
            int t = match[v];  
            if(t == 0 || find(t)) {  
                match[v] = u;  
                return true;  
            }  
        }  
    }  
    return false;  
}  	
for(int i = 1; i <= n; i++) {  
    memset(st, 0, sizeof st);  
    ans += find(i);  
}
```

### 最小路径覆盖 DAG

不相交

​	DAG，有向无环图，用最少的互不相交的路径，将所有点覆盖。

​	拆点，将原图中1...n，拆成新图1...n，1'...n'，原图i->j的有向边变成新图中i->j'的无向边，新图一定为二分图。

​	(实际上不需要新建一张图)

​	最小路径覆盖 = 原图点数 - 新图最大匹配数

可相交

​	用最少的可相交的路径，将所有点覆盖

​	对原图G，先求出传递闭包G'，然后对G’转化为不相交问题

code

```c
//floyd求传递闭包
for(int k = 1; k <= n; k ++)
        for(int i = 1; i <= n; i ++)
            for(int j = 1; j <= n; j ++)
                d[i][j] |= d[i][k] & d[k][j];
int ans = 0;
//匈牙利算法
for(int i = 1; i <= n; i ++) {
    memset(st, 0, sizeof st);
    ans += find(i);
}

cout << n - ans << "\n";
```

## 拓扑排序

### 无向图求环

```C
void topo()
{
	queue<int> q;
	
	for(int i=1;i<=n;i++)
		if(d[i]==1) q.push(i); //不是环内的点，入队
	
	while(!q.empty())
	{
		int u=q.front();
		vis[u]=1;	
		q.pop();
		for(int i=0;i<g[u].size();i++)
		{
			int v=g[u][i];
			if(--d[v]==1) q.push(v); //u和v断开，v的度数-1，减一后若度数变为，则可以入队
		}
	}
}
// 度数为2的点是环内的点
```

## 欧拉回路和欧拉路径

### 无向图

​	对所有边都连通的无向图

​	(1)存在欧拉路径的充分必要条件：度数为奇数的点只能有0或2个。

​	(2)存在欧拉回路的充分必要条件：度数为奇数的点只能有0个。

```c
void dfs(int u) {
    for(int i = 1; i <= 500; i ++) {
        if(g[u][i]) {
            g[u][i] --, g[i][u] --;
            dfs(i);
        }
    }
    ans[++ cnt] = u;
}
```

### 有向图

​	对所有边都连通的有向图

​	(1)存在欧拉路径的充分必要条件：要么所有点的出度均等于入度；要么除了两个点之外，其余所有点的出度等于入度，剩余的两个点：一个满足出度比入度多1（起点)，另一个满足入度比出度多1（终点）。

​	(2)存在欧拉回路的充分必要条件：所有点的出度均等于入度。

# 数据结构

## 并查集

```c
int get(int x) {
	if(fa[x] != x) fa[x] = get(fa[x]);
	return fa[x];
}

void merge(int x, int y) {
	x = get(x), y = get(y);
	if(x == y) return;
	fa[x] = y;
	val[y] += val[x];
}
```

## 树状数组

```c
int lowbit(int x)
{
    return x & -x;
}

void update(int x, int c)  // 位置x加c
{
    for (int i = x; i <= n; i += lowbit(i)) tr[i] += c;
}

int ask(int x) {
    int res = 0;
    for (int i = x; i; i -= lowbit(i)) res += tr[i];
    return res;
}
```

## 线段树

```c
struct Node
{
    int l, r;
    // TODO: 需要维护的信息和懒标记
}tr[N * 4];

void pushup(int u)
{
    // tr[u].v = max(tr[u << 1].v, tr[u << 1 | 1].v);
    // TODO: 利用左右儿子信息维护当前节点的信息
}

void pushdown(int u)
{
    // TODO: 将懒标记下传
}

void build(int u, int l, int r)
{
    if (l == r) tr[u] = {l, r};
    else
    {
        tr[u] = {l, r};
        int mid = l + r >> 1;
        build(u << 1, l, mid), build(u << 1 | 1, mid + 1, r);
        pushup(u);
    }
}

void update(int u, int l, int r, int d)
{
    if (tr[u].l >= l && tr[u].r <= r)
    {
        // TODO: 修改区间
    }
    else
    {
        pushdown(u);
        int mid = tr[u].l + tr[u].r >> 1;
        if (l <= mid) update(u << 1, l, r, d);
        if (r > mid) update(u << 1 | 1, l, r, d);
        pushup(u);
    }
}

int query(int u, int l, int r)
{
    if (tr[u].l >= l && tr[u].r <= r)
    {
        return ;  // TODO 需要补充返回值
    }
    else
    {
        pushdown(u);
        int mid = tr[u].l + tr[u].r >> 1;
        int res = 0;
        if (l <= mid ) res = query(u << 1, l, r);
        if (r > mid) res += query(u << 1 | 1, l, r);
        return res;
    }
}
```

### 扫描线

```c
void pushup(int u)
{
    if (tr[u].cnt) tr[u].len = ys[tr[u].r + 1] - ys[tr[u].l];
    else if (tr[u].l != tr[u].r)
    {
        tr[u].len = tr[u << 1].len + tr[u << 1 | 1].len;
    }
    else tr[u].len = 0;
}

void modify(int u, int l, int r, int k)
{
    if (tr[u].l >= l && tr[u].r <= r)
    {
        tr[u].cnt += k;
        pushup(u);
    }
    else
    {
        int mid = tr[u].l + tr[u].r >> 1;
        if (l <= mid) modify(u << 1, l, r, k);
        if (r > mid) modify(u << 1 | 1, l, r, k);
        pushup(u);
    }
}
int main() {
    for (int i = 0, j = 0; i < n; i ++ )
    {
        double x1, y1, x2, y2;
        scanf("%lf%lf%lf%lf", &x1, &y1, &x2, &y2);
        seg[j ++ ] = {x1, y1, y2, 1};
        seg[j ++ ] = {x2, y1, y2, -1};
        ys.push_back(y1), ys.push_back(y2);
    }

    sort(ys.begin(), ys.end());

    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    build(1, 0, ys.size() - 2);

    sort(seg, seg + n * 2);

    double res = 0;
    for (int i = 0; i < n * 2; i ++ )
    {
        if (i > 0) res += tr[1].len * (seg[i].x - seg[i - 1].x);
        modify(1, find(seg[i].y1), find(seg[i].y2) - 1, seg[i].k);
    }
}

```

## 分块

```c
int id[N], len, tot;
struct node {
	int l, r;
	ll sum, tag;
}blk[N];
int n, m;

void update(int l, int r, ll val) {
	int idx = id[l], idy = id[r];
	if(idx == idy) {
		for(int i = l; i <= r; i ++)
			a[i] += val;
		blk[idx].sum += (r - l + 1) * val;
	} else {
		for(int i = l; i <= blk[idx].r; i ++)
			a[i] += val;
		blk[idx].sum += (blk[idx].r - l + 1) * val;
		for(int i = blk[idy].l; i <= r; i ++) 
			a[i] += val;
		blk[idy].sum += (r - blk[idy].l + 1) * val;
		for(int i = idx + 1; i < idy;  i ++) {
			blk[i].tag += val;
			blk[i].sum += val * (blk[i].r - blk[i].l + 1);
		}
	}
}

ll query(int l, int r) {
	int idx = id[l], idy = id[r];
	ll ans = 0;
	if(idx == idy) {
		for(int i = l; i <= r; i ++) {
			ans += a[i] + blk[idx].tag;
		}
		return ans;
	} else {
		for(int i = l; i <= blk[idx].r; i ++)
			ans += a[i] + blk[idx].tag;
		for(int i = blk[idy].l; i <= r; i ++) 
			ans += a[i] + blk[idy].tag;
		for(int i = idx + 1; i < idy;  i ++)
			ans += blk[i].sum;
		return ans;
	}
}

void build() {
	len = sqrt(n);
	for(int i = 1; i <= n; i ++) {
		id[i] = (i - 1) / len + 1;
		blk[id[i]].sum += a[i];
		if((i - 1) % len==0) blk[++tot].l = i;
		if(i%len == 0) blk[tot].r = i;
	}
    blk[tot].r = n;
}
```

## 可持久化tire

```c++
int n, m;
int s[N];
int tr[M][2], max_id[M];
int root[N], idx;

void insert(int i, int k, int p, int q) {
    if(k < 0) {
        max_id[q] = i;
        return;
    }
    int v = s[i] >> k&1;
    if (p) tr[q][v ^ 1] = tr[p][v ^ 1];
    tr[q][v] = ++ idx;
    insert(i, k - 1, tr[p][v], tr[q][v]);
    max_id[q] = max(max_id[tr[q][0]], max_id[tr[q][1]]);
}

int query(int root, int C, int L) {
    int p = root;
    for (int i = 23; i >= 0; i --) {
        int v = C >> i & 1;
        if (max_id[tr[p][v ^ 1]] >= L) p = tr[p][v ^ 1];
        else p = tr[p][v];
    }
    return C ^ s[max_id[p]];
}

void tire_init() {
    max_id[0] = -1;
    root[0] = ++idx;
    insert(0, 23, 0, root[0]);
}

for(int i = 1; i <= n; i ++) {
        ll tmp;
        cin >> tmp;
        s[i] = s[i - 1] ^ tmp;
        root[i] = ++ idx;
        insert(i, 23, root[i - 1], root[i]);
   }
```



# 字符串

## AC自动机

###  一

```c
//有多少个不同的模式串在文本串里出现过 
namespace AC {
    int tr[N][26], tot;
    int e[N], fail[N];

    void insert(char *s) {
        int u = 0;
        for (int i = 1; s[i]; i ++) {
            if(!tr[u][s[i] - 'a']) tr[u][s[i] - 'a'] = ++tot;
            u = tr[u][s[i] - 'a'];
        }
        e[u] ++;
    }

    queue<int> q;

    void build() {
        for(int i = 0; i < 26; i++)
            if(tr[0][i]) q.push(tr[0][i]);
        while (q.size()) {
            int u = q.front();
            q.pop();
            for (int i = 0; i < 26; i ++) {
                if (tr[u][i]) {
                    fail[tr[u][i]] = tr[fail[u]][i];
                    q.push(tr[u][i]);
                } else tr[u][i] = tr[fail[u]][i];
            }
        } 
    }

    int query(char *t) {
        int u = 0, res = 0;
        for (int i = 1; t[i]; i++) {
            u = tr[u][t[i] - 'a'];
            for (int j = u; j && e[j] != -1; j = fail[j]) {
                res += e[j], e[j] = -1;
            }
        }
        return res;
    }
}
```

### 二

```c
//AC自动机
#include <bits/stdc++.h>
using namespace std;

const int N = 2e6 + 10;
char s[N];
char a[N];
namespace AC {
    int tr[N][26], tot;
    int idx[N], fail[N];
    int val[N];
    int vis[N];
    int cnt[N]; // 第i个字符串的出现次数

    void insert(char *s, int id) {
        int u = 0;
        for (int i = 1; s[i]; i ++) {
            if(!tr[u][s[i] - 'a']) tr[u][s[i] - 'a'] = ++tot;
            u = tr[u][s[i] - 'a'];
        }
        vis[id] = u;
    }

    void init() {
        memset(fail, 0, sizeof fail);
        memset(tr, 0, sizeof tr);
        memset(idx, 0, sizeof idx);
        memset(val, 0, sizeof val);
        memset(cnt, 0, sizeof cnt);
        tot = 0;
    }

    queue<int> q;

    void build() {
        for(int i = 0; i < 26; i++)
            if(tr[0][i]) q.push(tr[0][i]);
        while (q.size()) {
            int u = q.front();
            q.pop();
            for (int i = 0; i < 26; i ++) {
                if (tr[u][i]) {
                    fail[tr[u][i]] = tr[fail[u]][i];
                    q.push(tr[u][i]);
                } else tr[u][i] = tr[fail[u]][i];
            }
        } 
    }

    int query(char *t) {
        int u = 0, res = 0;
        for (int i = 1; t[i]; i++) {
            u = tr[u][t[i] - 'a'];
            for (int j = u; j ; j = fail[j]) res val[j] ++;
        }
        return res;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    int n;
    cin >> n;
    // AC::init();
    for(int i = 1; i <= n; i ++) {
        cin >> s + 1;
        AC::insert(s, i);
    }
    AC::build();
    cin >> a + 1;
    int x = AC::query(a);
    for(int i = 1; i <= n; i ++) {
        cout << AC::val[AC::vis[i]] << "\n";
    }
}
```

### 三.树上差分优化

```c
    int head[N], nxt[N], to[N], cnt;

    void add(int u, int v) {
            nxt[++cnt] = head[u];
            head[u] = cnt;
            to[cnt] = v;
    }
    
    void dfs(int u) {
        int i, v;
        for(i = head[u]; i; i = nxt[i]) {
            v = to[i];
            dfs(v);
            val[u] += val[v];
        }
    }

    int query(char *t) {
        int u = 0, res = 0;
        for (int i = 1; t[i]; i++) {
            u = tr[u][t[i] - 'a'];
            ++val[u];
        }

        for(int i = 1; i <= tot; i ++) add(fail[i], i);

        dfs(0);

        return res;
    }

```

## 后缀自动机

### 节点含义

​	每个节点代表的是一个连续长度的字符串的集合，这个集合里面的字符串的长度是连续的，并且每个长度都只会出现一次，字符串集合里面字符串都是长度比它大的字符串的后缀

### Fail指针

​	跳到比当前节点短的最长的相同后缀字符串集

​	Fail[x]的字符串集长度和x的字符串集长度是连续的

```c
struct Suffix_Automata {
    //两倍字符串长度的空间
    int maxlen[Maxn], trans[Maxn][26], link[Maxn], Size, Last;
    int a[Maxn], b[Maxn], endpos[Maxn];
    int dp[Maxn], sum[Maxn];
    Suffix_Automata() {Size = Last = 1;}

    inline void Extend(int id) {
        int cur = (++ Size), p;
        maxlen[cur] = maxlen[Last] + 1;
        endpos[cur] = 1;
        for(p = Last; p && !trans[p][id]; p = link[p]) trans[p][id] = cur;
        if (!p) link[cur] = 1;
        else {
            int q = trans[p][id];
            if (maxlen[q] == maxlen[p] + 1) link[cur] = q;
            else {
                int clone = (++ Size);
                maxlen[clone] = maxlen[p] + 1;
                for(int i = 0; i < 26; ++i) trans[clone][i] = trans[q][i];
                link[clone] = link[q];
                for(; p && trans[p][id] == q; p = link[p]) trans[p][id] = clone;
                link[cur] = link[q] = clone;
            }
        }
        Last = cur;
    }

    void getTP(int &Len) { //getendpos 前置
        for(int i = 1; i <= Size; i ++) a[maxlen[i]] ++;
        for(int i = 1; i <= Len; i ++) a[i] += a[i - 1];
        for(int i = 1; i <= Size; i++) b[a[maxlen[i]] --] = i;
    }

    void getendpos() { //求每类子串的数量
        for(int i = Size; i >= 1; i --) {
            int e = b[i];
            endpos[link[e]] += endpos[e];
        }
    }

    ll getSubNum() { // 求不相同子串数量
        ll ans = 0;
        for(int i = 2 ;i <= Size; ++i)
            ans += maxlen[i] - maxlen[link[i]] ;
        return ans;
    }

    void get_len_max(int Len) { // 求长度为i的出现次数最多的子串
        for(int i = 1; i <= Size; i++) dp[maxlen[i]] = max(dp[maxlen[i]], endpos[i]);
        for(int i = Len - 1; i >= 1; i --) dp[i] = max(dp[i], dp[i + 1]);
        for(int i = 1; i <= Len; i ++) {
            cout << dp[i] << "\n";
        }
    }

    void getsum() { //getmink前置
        for(int i = Size; i >= 1; i --) {
            int &e = b[i];
            sum[e] = 1;
            for(int j = 0; j < 26; j ++) {
                sum[e] += sum[trans[e][j]];
            }
        }
    }

    void getmink(int k) { // 求字典序第k小的子串
        int now = 1, p;
        string s = "";
        while(k) {
            for (int i = 0; i < 26; i ++) {
                if (trans[now][i] && k) {
                    p = trans[now][i];
                    if (sum[p] < k) k -= sum[p];
                    else {
                        k --; now = p;
                        s += (char)(i + 'a');
                        break;
                    }
                }
            }
        }
        cout << s << "\n";
    }

    void LCS(char s[], int Len) { // 求两个串的最长公共子串
        int ans = 0, cnt = 0;
        int now = 1;
        char base = 'a';
        for(int i = 0; i < Len; i++) {
            int c = s[i] - base;
            if (trans[now][c]) {
                cnt ++;
                now = trans[now][c];
            } else {
                while(now&&!trans[now][c]) now = link[now];
                if (!now) cnt=0, now=1;
                else cnt = maxlen[now]+1, now = trans[now][c];
            }
            ans = max(ans, cnt);
        }
        cout << ans << "\n";
    }

    int ans[Maxn];
    void init_ans(){
        for(int i=1;i<=Size;i++) ans[i]= maxlen[i];
    }
    
    void LCS2(char s[],int Len){     //求多个串的最长公共子串 
        for(int i=1;i<=Size;i++) dp[i]=0;
        int cnt=0;
        int now=1;
        char base='a';
        for(int i=0;i<Len;i++){
            int c=s[i]-base;
            if(trans[now][c]){
                cnt++;
                now=trans[now][c];
            }
            else{
                while(now&&!trans[now][c]) now=link[now];
                if(!now) cnt=0,now=1;
                else cnt=maxlen[now]+1,now=trans[now][c];
            }
            dp[now]=max(dp[now],cnt);
        }
        for(int i=Size;i>=1;i--){
            int e=b[i];
            dp[link[e]]=max(dp[link[e]],min(dp[e],maxlen[link[e]]));
        }
        for(int i=1;i<=Size;i++) ans[i]=min(ans[i],dp[i]);
    }
    void get_LCS2_ans(){
        int cnt=0;
        for(int i=1;i<=Size;i++) cnt=max(cnt,ans[i]);
        printf("%d\n",cnt);
    }

    void get_cntk(int k) { //求出现次数为k的子串种数
        ll ans = 0;
        for (int i = 1; i <= Size; i++) {
            if(endpos[i] == k) ans += maxlen[i] - maxlen[link[i]];
        }
        cout << ans << "\n";
    }

    int d[Maxn];
    void get_sumk(int l, int r) { //求出现次数 l <= k <= r 的子串种数
        for(int i = Size; i > 1; i--) {
            int v = b[i];
            if(endpos[v] >= l && endpos[v] <= r) d[v] ++;
            for(int j = 0; j < 26; j ++) {
                if(trans[v][j]) d[v] += d[trans[v][j]];
            }
        }
        ll ans = 0;
        for (int i = 0; i < 26; i ++) {
            if (trans[1][i]) ans += d[trans[1][i]];
        }
        cout << ans << "\n";
    }
} T;
```

# 组合数学

## 组合

### 可重集合的组合公式

​	有k种元素，每个元素有ai，从中无序地选出r个元素，r <= ai，则有选法：
$$
C_{k-1+r}^{r}
$$


## 排列

### 可重集合的排列公式

​	有k个元素，每个元素有ai个，不同排列的个数为
$$
\frac{(\sum_{i=1}^{k}a_i)!}{\prod_{i=1}^k(a_i!)}
$$

# 常见技巧

## python

### 输入与输出

第一行包含一个整数 n，表示待处理的整数个数。

第二行包含空格分隔的 n 个整数，依次表示 a1,a2,⋯,an。

```python
n = int(input())
arr = [int(x) for x in sys.stdin.readline().split()]
```

### map

```python
# 初始化
vis = {}
# 是否在vis内
vis.has_key(key)
# 以列表返回可遍历的(键, 值) 元组数组
vis.items()
# 以列表返回一个字典所有的键
vis.keys()
# 删除字典内所有元素
vis.clear()
# 删除字典给定键key对应的值 
vis.pop(key)
```

