//=======================================
// Chapter 4, RBC model, deterministic solution
//   modified on 2020/01/02
//=======================================

// 1. 内生変数、外生変数の宣言
var C L K Y w r A;
varexo e;

// 2. パラメータの宣言
parameters alpha beta delta mu gamma rho;

// パラメータ値の代入
alpha = 0.3;
beta = 0.99;
delta = 0.025;
mu = 1.0;
gamma = 1.0;
rho = 0.9;

// 3. 方程式の定義
model;
w/C = (gamma+1)*mu*L^gamma;
C(+1)/C = beta*(r(+1)-delta+1);
Y = A*K^alpha*L^(1-alpha);
w = (1-alpha)*A*K^alpha*L^(-alpha);
r = alpha*A*K^(alpha-1)*L^(1-alpha);
K = Y(-1)+(1-delta)*K(-1)-C(-1);
log(A) = rho*log(A(-1)) + e;
end;

// 4. 定常状態の計算
Astar = 1;
rstar = 1/beta + delta - 1;
K_L = (rstar/alpha/Astar)^(1/(alpha-1));
Y_L = Astar*K_L^alpha;
C_L = Y_L-delta*K_L;
wstar = (1-alpha)*Astar*K_L^alpha;
Lstar = (wstar/(gamma+1)/mu)^(1/(gamma+1))*C_L^(-1/(gamma+1));
Kstar = K_L*Lstar;
Ystar = Y_L*Lstar;
Cstar = C_L*Lstar;

// Dynareに定常状態を計算させる際の初期値
initval;
C = Cstar;
L = Lstar;
K = Kstar;
Y = Ystar;
w = wstar;
r = rstar;
A = Astar;
end;

// Dynareに定常状態を計算させる
steady;

// 初期値と結果が変わらないことをチェック
[Cstar; Lstar; Kstar; Ystar; wstar; rstar; Astar]

// モデルのチェック
check;

// 5. シミュレーション
// シナリオの設定
shocks;
var e;
periods 1;
values 0.01;
end;

// シミュレーションの実行
simul(periods=150);

// 定常状態からの乖離率の計算
C1 = (C./Cstar-1)*100;
L1 = (L./Lstar-1)*100;
K1 = (K./Kstar-1)*100;
Y1 = (Y./Ystar-1)*100;
w1 = (w./wstar-1)*100;
r1 = (r-rstar)*100;
A1 = (A./Astar-1)*100;
I1 = ((Y-C)./(Ystar-Cstar)-1)*100;

// グラフ描写
figure(1)
subplot(2,2,1)
plot(0:50, A1(1:51)); title('A')
subplot(2,2,2)
plot(0:50, Y1(1:51)); title('Y')
subplot(2,2,3)
plot(0:50, C1(1:51)); title('C')
subplot(2,2,4)
plot(0:50, K1(1:51)); title('K')

figure(2)
subplot(2,2,1)
plot(0:50, L1(1:51)); title('L')
subplot(2,2,2)
plot(0:50, I1(1:51)); title('I')
subplot(2,2,3)
plot(0:50, w1(1:51)); title('w')
subplot(2,2,4)
plot(0:50, r1(1:51)); title('r')
