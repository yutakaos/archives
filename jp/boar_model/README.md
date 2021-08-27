# 解説記事

このページでは和文記事の JAGS コードを説明しています。
論文は J-STAGE からダウンロードができます。
https://www.jstage.jst.go.jp/article/hozen/23/1/23_29/_article/-char/ja/

＊説明するにあたり、出版当時に使用したコードとは同じでなく改変しています。

## データ構造（data）の指定

値を推定しない以下の変数を指定します。

1. 変数選択の事前分布
2. One's trick と Zero's trick のためのダミー変数

```r
data {
    # 1. Uninformative priors for variable selection
    for (k in 1:NV) {
        PI[k] <- 0.5
    }
    
    # 2. Dummy variables for one's and zero's trick
    for (i in 1:NS) {
        ONE[i,1] <- 1
        ONE[i,2] <- 1
        ONE[i,3] <- 1
        ZERO[i]  <- 0
    }
}
```

## One's trick と Zero's trick

本コードでは One's trick と Zero's trick を使用するため、簡単に説明します。

JAGS では、多くの主要な確率分布が事前に用意されており、
多くの場合に確率分布をチルダで指定すれば簡単に尤度を計算することができます。

1. 通常の方法

```r
model{
    Y ~ dbern(p)
}
```

しかしながら、事前に用意されていない確率分布には通常の方法は使えません。
One's trick や Zero's trick は通常とは異なる記法で尤度を計算するためのテクニックです。

2. One's trick（L = ベルヌーイ分布の尤度）

```r
data {
    ONE <- 1
}
model {
    L <- p^Y + (1-p)^(1-Y)
    ONE ~ dbern(L)
}
```

3. Zero's trick（logL = ベルヌーイ分布の対数尤度）

```r
data {
    ZERO <- 0
}
model {
    logL <- Y * log(p) + (1-Y) * log(1-p)
    ZERO ~ dpois(-logL)
}
```

確率分布によって尤度が計算しやすい場合と対数尤度が計算しやすい場合があるので、
状況によって One's trick と Zero's trick を使い分けます。


## モデル構造（model）の指定

モデルを書くにあたり、「データなどの値が一定の変数」は大文字、「推定パラメータなど値の変化する変数」は 小文字というルールで書いています。


### 観測モデル

- 観測モデルでは、市町単位の箱わな捕獲数（CAPT）がポアソン分布に従うことを仮定しています。

平均捕獲数（mu） = 捕獲率（r_capt）× 生息個体数（num_g）/ 生息地面積（AREA）× 箱わな数（TRAP）

- 市町・年によっては捕獲数あるいは箱わな数のデータが書けているため、
NCS と IDCS という変数を使って欠損値を無視しています。

```r
# 1. Observation model
for (y in 1:2) {
    for (i in 1:NCS[y]) {
        mu[i,y] <- r_capt * num_g[IDCS[i,y],y] / AREA[IDCS[i,y]] * TRAP[IDCS[i,y],y]
        CAPT[IDCS[i,y],y] ~ dpois(mu[i,y])
    }
}
r_capt ~ dunif(0, 5)
```

### システムモデル

- システムモデルでは、市町単位のイノシシ個体群が

　　・・・ → 捕獲による個体数減少 → 出生死亡による個体数増加 → ・・・

というプロセスにより個体群変動することを想定しています。

- 実際には移入出も個体群変動に影響しているはずですが、
市町単位での解析は移入出より出生死亡の影響が大きいと考えられます。
また、個体群成長率（growth）には変量効果が含まれているので、
移入出の影響があれば変量効果によって考慮されます。

- 捕獲前の個体数（num_g）は捕獲数（HUNT）より小さくなりません。
そこで、One's trick により HUNT < num_g であることを保証します。

```r
# 2. System model
for (i in 1:NS) {
    # Convert to density
    density[i,1] <- num_h[i,1] / AREA[i]
    density[i,2] <- num_h[i,2] / AREA[i]
    
    # Population dynamics
    # y = 1
    num_0[i] ~ dunif(1,10000)
    num_h[i,1] <- trunc(num_0[i])
    num_g[i,1] ~ dpois(growth[i] * num_h[i,1])
    
    # y = 2
    num_h[i,2] <- num_g[i,1] - HUNT[i,1]
    num_g[i,2] ~ dpois(growth[i] * num_h[i,2])
    
    # Constraints (HUNT < num_g)
    ONE[i,1] ~ dbern(step(num_g[i,1] - HUNT[i,1]))
    ONE[i,2] ~ dbern(step(num_g[i,2] - HUNT[i,2]))
}
```

### パラメータモデル

- 個体群成長率を線形混合モデル（LMM）によってモデル化します。

- 個体群成長率は 8.0 以下になるよう制約しています。

- 変量効果の事前分布は半コーシー分布を使っています。

```r
# 3. Parameter model
# Priors
g_prec ~ dgamma(0.001, 0.001)
prec[1] <- g_prec
prec[2] <- g_prec * (sum(vs) + 1)

# Intercept and slops (with variable selection)
g_coef[1] ~ dnorm(0, prec[2])
for (k in 1:NV) {
    vs  [k] ~ dbern(PI[k])
    coef[k] ~ dnorm(0, prec[vs[k] + 1])
    g_coef[k+1] <- coef[k] * vs[k]
}

# Priors for random effects (half-cauchy)
tau <- 1 / sigma / sigma
sigma ~ dt(0, 0.0016, 1) I(0,)

# Growth rate
for (i in 1:NS) {
    # Linear mixed model (LMM)
    log(growth[i]) <- inprod(g_coef, ENVI[i,]) + g_rand[i]
    g_rand[i] ~ dnorm(0.0, tau)
    
    # Constraints (growth < 8.0)
    ONE[i,3] ~ dbern(step(8.0 - growth[i]))
}
```

### 初年度個体数の空間自己相関

- 調査地千葉県のイノシシは歴史的な分布拡大をしてきた経緯があります。
そこで、歴史的な密度の類似性を表現するため、初年度に空間自己相関を考えます。

- s_mu =（隣り合う市町の密度の平均）ー（その市町の密度）

- 印西市はとても最近にイノシシの人為的侵入があったと考えられており、
千葉県北部で唯一イノシシが定着している地域でした（調査当時）。
そのため、印西市の空間自己相関は NSS と IDSS という変数を使って無視しています。

```r
# 4. spatial autocorrelation

for (i in 1:NSS) {
    s_mu[i] <- inprod(density[,1], WEIS[,IDSS[i]])
    ZERO[i] ~ dnorm(s_mu[i], s_prec)
}
s_prec ~ dgamma(0.001, 0.001)
```

- Zero's trick によっても同様に表現できます。

```r
# 4. spatial autocorrelation
for (i in 1:NSS) {
    s_mu[i] <- inprod(density[,1], WEIS[,IDSS[i]])
    zero[i] <- 0.92 + 0.5 * (-log(s_prec) + pow(s_mu[i], 2) * s_prec)
    ZERO[i] ~ dpois(zero[i])
}
s_prec ~ dgamma(0.001, 0.001)
```
