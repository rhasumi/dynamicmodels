//=======================================
// Chapter 5, New Keynesian model,
//      linearized, deterministic solution
//   modified on 2020/01/02
//=======================================

// 1. 内生変数、外生変数の宣言
var x ppi a ii v;
varexo e z;

// 2. パラメータの宣言
parameters beta kappa phi_pi phi_y rho_A rho_v;

// パラメータ値の代入
beta = 0.99;
gamma = 5;
varrho = 0.9;
kappa = (1-varrho)*(1-varrho*beta)*(gamma+1)/varrho;
phi_pi = 1.5;
phi_y = 0.5;
rho_A = 0.9;
rho_v = 0.7;

// 3. 方程式の定義
model(linear);
ppi = beta*ppi(+1)+kappa*x;
x = x(+1)-(ii-ppi(+1))+(rho_A-1)*a;
ii = phi_pi*ppi + phi_y*x + v;
a = rho_A*a(-1) + e;
v = rho_v*v(-1) + z;
end;

// 全てゼロになることを念のためチェック
steady;

// モデルのチェック
check;

// 5. シミュレーション(deterministic)
shocks;
var e; periods 1; values 0.01;
end;

simul(periods=150);

figure(1)
subplot(3,2,1)
plot(0:50, x(1:51)*100); title('x')
subplot(3,2,2)
plot(0:50, ii(1:51)*100); title('i')
subplot(3,2,3)
plot(0:50, ppi(1:51)*100); title('pi')
subplot(3,2,5)
plot(0:50, a(1:51)*100); title('a')
subplot(3,2,6)
plot(0:50, v(1:51)*100); title('v')

shocks;
var e; periods 1; values 0;
var z; periods 1; values 0.01;
end;

simul(periods=150);

figure(2)
subplot(3,2,1)
plot(0:50, x(1:51)*100); title('x')
subplot(3,2,2)
plot(0:50, ii(1:51)*100); title('i')
subplot(3,2,3)
plot(0:50, ppi(1:51)*100); title('pi')
subplot(3,2,5)
plot(0:50, a(1:51)*100); title('a')
subplot(3,2,6)
plot(0:50, v(1:51)*100); title('v')
