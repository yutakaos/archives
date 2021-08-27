# ����L��

���̃y�[�W�ł͘a���L���� JAGS �R�[�h��������Ă��܂��B
�_���� J-STAGE ����_�E�����[�h���ł��܂��B
https://www.jstage.jst.go.jp/article/hozen/23/1/23_29/_article/-char/ja/

����������ɂ�����A�o�œ����Ɏg�p�����R�[�h�Ƃ͓����łȂ����ς��Ă��܂��B

## �f�[�^�\���idata�j�̎w��

�l�𐄒肵�Ȃ��ȉ��̕ϐ����w�肵�܂��B

1. �ϐ��I���̎��O���z
2. One's trick �� Zero's trick �̂��߂̃_�~�[�ϐ�

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

## One's trick �� Zero's trick

�{�R�[�h�ł� One's trick �� Zero's trick ���g�p���邽�߁A�ȒP�ɐ������܂��B

JAGS �ł́A�����̎�v�Ȋm�����z�����O�ɗp�ӂ���Ă���A
�����̏ꍇ�Ɋm�����z���`���_�Ŏw�肷��ΊȒP�ɖޓx���v�Z���邱�Ƃ��ł��܂��B

1. �ʏ�̕��@

```r
model{
    Y ~ dbern(p)
}
```

�������Ȃ���A���O�ɗp�ӂ���Ă��Ȃ��m�����z�ɂ͒ʏ�̕��@�͎g���܂���B
One's trick �� Zero's trick �͒ʏ�Ƃ͈قȂ�L�@�Ŗޓx���v�Z���邽�߂̃e�N�j�b�N�ł��B

2. One's trick�iL = �x���k�[�C���z�̖ޓx�j

```r
data {
    ONE <- 1
}
model {
    L <- p^Y + (1-p)^(1-Y)
    ONE ~ dbern(L)
}
```

3. Zero's trick�ilogL = �x���k�[�C���z�̑ΐ��ޓx�j

```r
data {
    ZERO <- 0
}
model {
    logL <- Y * log(p) + (1-Y) * log(1-p)
    ZERO ~ dpois(-logL)
}
```

�m�����z�ɂ���Ėޓx���v�Z���₷���ꍇ�Ƒΐ��ޓx���v�Z���₷���ꍇ������̂ŁA
�󋵂ɂ���� One's trick �� Zero's trick ���g�������܂��B


## ���f���\���imodel�j�̎w��

���f���������ɂ�����A�u�f�[�^�Ȃǂ̒l�����̕ϐ��v�͑啶���A�u����p�����[�^�Ȃǒl�̕ω�����ϐ��v�� �������Ƃ������[���ŏ����Ă��܂��B


### �ϑ����f��

- �ϑ����f���ł́A�s���P�ʂ̔���ȕߊl���iCAPT�j���|�A�\�����z�ɏ]�����Ƃ����肵�Ă��܂��B

���ϕߊl���imu�j = �ߊl���ir_capt�j�~ �����̐��inum_g�j/ �����n�ʐρiAREA�j�~ ����Ȑ��iTRAP�j

- �s���E�N�ɂ���Ă͕ߊl�����邢�͔���Ȑ��̃f�[�^�������Ă��邽�߁A
NCS �� IDCS �Ƃ����ϐ����g���Č����l�𖳎����Ă��܂��B

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

### �V�X�e�����f��

- �V�X�e�����f���ł́A�s���P�ʂ̃C�m�V�V�̌Q��

�@�@�E�E�E �� �ߊl�ɂ��̐����� �� �o�����S�ɂ��̐����� �� �E�E�E

�Ƃ����v���Z�X�ɂ��̌Q�ϓ����邱�Ƃ�z�肵�Ă��܂��B

- ���ۂɂ͈ړ��o���̌Q�ϓ��ɉe�����Ă���͂��ł����A
�s���P�ʂł̉�͈͂ړ��o���o�����S�̉e�����傫���ƍl�����܂��B
�܂��A�̌Q�������igrowth�j�ɂ͕ϗʌ��ʂ��܂܂�Ă���̂ŁA
�ړ��o�̉e��������Εϗʌ��ʂɂ���čl������܂��B

- �ߊl�O�̌̐��inum_g�j�͕ߊl���iHUNT�j��菬�����Ȃ�܂���B
�����ŁAOne's trick �ɂ�� HUNT < num_g �ł��邱�Ƃ�ۏ؂��܂��B

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

### �p�����[�^���f��

- �̌Q����������`�������f���iLMM�j�ɂ���ă��f�������܂��B

- �̌Q�������� 8.0 �ȉ��ɂȂ�悤���񂵂Ă��܂��B

- �ϗʌ��ʂ̎��O���z�͔��R�[�V�[���z���g���Ă��܂��B

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

### ���N�x�̐��̋�Ԏ��ȑ���

- �����n��t���̃C�m�V�V�͗��j�I�ȕ��z�g������Ă����o�܂�����܂��B
�����ŁA���j�I�Ȗ��x�̗ގ�����\�����邽�߁A���N�x�ɋ�Ԏ��ȑ��ւ��l���܂��B

- s_mu =�i�ׂ荇���s���̖��x�̕��ρj�[�i���̎s���̖��x�j

- �󐼎s�͂ƂĂ��ŋ߂ɃC�m�V�V�̐l�דI�N�����������ƍl�����Ă���A
��t���k���ŗB��C�m�V�V���蒅���Ă���n��ł����i���������j�B
���̂��߁A�󐼎s�̋�Ԏ��ȑ��ւ� NSS �� IDSS �Ƃ����ϐ����g���Ė������Ă��܂��B

```r
# 4. spatial autocorrelation

for (i in 1:NSS) {
    s_mu[i] <- inprod(density[,1], WEIS[,IDSS[i]])
    ZERO[i] ~ dnorm(s_mu[i], s_prec)
}
s_prec ~ dgamma(0.001, 0.001)
```

- Zero's trick �ɂ���Ă����l�ɕ\���ł��܂��B

```r
# 4. spatial autocorrelation
for (i in 1:NSS) {
    s_mu[i] <- inprod(density[,1], WEIS[,IDSS[i]])
    zero[i] <- 0.92 + 0.5 * (-log(s_prec) + pow(s_mu[i], 2) * s_prec)
    ZERO[i] ~ dpois(zero[i])
}
s_prec ~ dgamma(0.001, 0.001)
```
