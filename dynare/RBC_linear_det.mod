//=======================================
// Chapter 4, RBC model, 
//     linearized, deterministic solution
//   modified on 2021/03/03
//=======================================

// 1. 内生変数、外生変数の宣言
var C L K Y w R A;
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
w-C = gamma*L;
C(+1)-C = beta*Rstar*R(+1);
Y = A + alpha*K + (1-alpha)*L;
w = A + alpha*K - alpha*L;
Rstar/(Rstar-1)*R = A + (alpha-1)*K + (1-alpha)*L;
Kstar*K = Ystar*Y(-1)+(1-delta)*Kstar*K(-1)-Cstar*C(-1);
A = rho*A(-1) + e;
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
//simul(periods=150);
perfect_foresight_setup(periods=150);
perfect_foresight_solver;
for i = 1:size(M_.endo_names, 1)
  assignin('base', string(M_.endo_names(i)), oo_.endo_simul(i,:)');
end

I = (Y*Ystar-C*Cstar)./(Ystar-Cstar);

// グラフ描写
figure(3)
subplot(2,2,1)
plot(0:50, A(1:51)*100); title('A')
subplot(2,2,2)
plot(0:50, Y(1:51)*100); title('Y')
subplot(2,2,3)
plot(0:50, C(1:51)*100); title('C')
subplot(2,2,4)
plot(0:50, K(1:51)*100); title('K')

figure(4)
subplot(2,2,1)
plot(0:50, L(1:51)*100); title('L')
subplot(2,2,2)
plot(0:50, I(1:51)*100); title('I')
subplot(2,2,3)
plot(0:50, w(1:51)*100); title('w')
subplot(2,2,4)
plot(0:50, R(1:51)*100); title('r')
