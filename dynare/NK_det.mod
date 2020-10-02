//=======================================
// Chapter 5, New Keynesian model,
//      deterministic solution
//   modified on 2020/01/02
//=======================================

// 1. 内生変数、外生変数の宣言
var C phi ppi ppitil F D A ii v;
varexo e z;

// 2. パラメータの宣言
parameters beta mu gamma varrho eta phi_pi phi_y rho_A rho_v istar Cstar;

// パラメータ値の代入
beta = 0.99;
mu = 1.0;
gamma = 5;
eta = 10;
varrho = 0.9;
phi_pi = 1.5;
phi_y = 0.5;
rho_A = 0.9;
rho_v = 0.7;

// 3. 方程式の定義
model;
(1+ppi(+1))*C(+1)/C = beta*(1+ii);
(1+ppitil)/(1+ppi) = eta/(eta-1)*F/D;
F = phi + varrho*beta*(1+ppi(+1))^eta*F(+1);
D = 1 + varrho*beta*(1+ppi(+1))^(eta-1)*D(+1);
(1+ppi)^(1-eta) = (1-varrho)*(1+ppitil)^(1-eta)+varrho;
ii = istar + phi_pi*ppi + phi_y*log(C/Cstar/A) + v;
log(A) = rho_A*log(A(-1)) + e;
v = rho_v*v(-1) + z;
phi - mu*(C/A)^(gamma+1)*(gamma+1) = 0;
end;

// 4. 定常状態の計算
pistar = 0;
pitilstar = 0;
vstar = 0;
istar = (1-beta)/beta;
Dstar = 1/(1-varrho*beta);
Fstar = (eta-1)/eta/(1-varrho*beta);
phistar = (eta-1)/eta;
Astar = 1;
Cstar = ((eta-1)/eta/mu/(gamma+1))^(1/(gamma+1));

initval;
C  = Cstar;
phi = phistar;
ppi = pistar;
ppitil = pitilstar;
F = Fstar;
D = Dstar;
A = Astar;
ii  = istar;
v = vstar;
end;

steady;

// モデルのチェック
check;

// 5. シミュレーション(deterministic)
// シナリオの設定
shocks;
var e; periods 1; values 0.01;
end;

// シミュレーションの実行
simul(periods=150);

// グラフ描写など
C1 = (C./Cstar-1)*100;
i1 = (ii-istar)*100;
pi1 = (ppi-pistar)*100;
A1 = (A./Astar-1)*100;
phi1 = (phi./phistar-1)*100;
v1 = (v-vstar)*100;

figure(1)
subplot(3,2,1)
plot(0:50, C1(1:51)); title('C')
subplot(3,2,2)
plot(0:50, i1(1:51)); title('i')
subplot(3,2,3)
plot(0:50, pi1(1:51)); title('pi')
subplot(3,2,4)
plot(0:50, phi1(1:51)); title('phi')
subplot(3,2,5)
plot(0:50, A1(1:51)); title('A')
subplot(3,2,6)
plot(0:50, v1(1:51)); title('v')

shocks;
var e; periods 1; values 0;
var z; periods 1; values 0.01;
end;

simul(periods=150);

C2 = (C./Cstar-1)*100;
i2 = (ii-istar)*100;
pi2 = (ppi-pistar)*100;
A2 = (A./Astar-1)*100;
phi2 = (phi./phistar-1)*100;
v2 = (v-vstar)*100;

figure(2)
subplot(3,2,1)
plot(0:50, C2(1:51)); title('C')
subplot(3,2,2)
plot(0:50, i2(1:51)); title('i')
subplot(3,2,3)
plot(0:50, pi2(1:51)); title('pi')
subplot(3,2,4)
plot(0:50, phi2(1:51)); title('phi')
subplot(3,2,5)
plot(0:50, A2(1:51)); title('A')
subplot(3,2,6)
plot(0:50, v2(1:51)); title('v')
