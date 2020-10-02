//=======================================
// Chapter 4, RBC model, 
//     linearized, stochastic solution
//   modified on 2020/01/02
//=======================================

// 1. 内生変数、外生変数の宣言
var c l k y w R a;
varexo e;

// 2. パラメータの宣言
parameters alpha beta delta mu gamma rho Rstar Kstar Ystar Cstar;

// パラメータ値の代入
alpha = 0.3;
beta = 0.99;
delta = 0.025;
mu = 1.0;
gamma = 1.0;
rho = 0.9;

model;
w-c = gamma*l;
c(+1)-c = beta*Rstar*R(+1);
y = a + alpha*k + (1-alpha)*l;
w = a + alpha*k - alpha*l;
Rstar/(Rstar-1)*R = a + (alpha-1)*k + (1-alpha)*l;
Kstar*k = Ystar*y(-1)+(1-delta)*Kstar*k(-1)-Cstar*c(-1);
a = rho*a(-1) + e;
end;

// 4. 定常状態の計算
Astar = 1;
rstar = 1/beta+delta-1;
K_L = (rstar/alpha/Astar)^(1/(alpha-1));
Y_L = Astar*K_L^alpha;
C_L = Y_L-delta*K_L;
wstar = (1-alpha)*Astar*K_L^alpha;
Lstar = (wstar/(gamma+1)/mu)^(1/(gamma+1))*C_L^(-1/(gamma+1));
Kstar = K_L*Lstar;
Ystar = Y_L*Lstar;
Cstar = C_L*Lstar;
Rstar = rstar+1;

// 全てゼロになることを念のためチェック
steady;

// モデルのチェック
check;

// 5. シミュレーション(stochastic)
// シナリオの設定
shocks;
var e = 1;
end;

// シミュレーションの実行
stoch_simul(order=1, irf = 100) c l k y w R a;
