//=======================================
// Chapter 4, RBC model, 
//     linearized, deterministic solution
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

figure(2)
subplot(2,2,1)
plot(0:50, L1(1:51)); title('L')
subplot(2,2,2)
plot(0:50, I1(1:51)); title('I')
subplot(2,2,3)
plot(0:50, w1(1:51)); title('w')
subplot(2,2,4)
plot(0:50, r1(1:51)); title('r')
