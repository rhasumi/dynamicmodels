//===========================================
// Chapter 8, Model definition and estimation
//   modified on 2019/10/19
//===========================================

// 1. 内生変数、外生変数の宣言
var x ppi a ii v;
varexo epsilon z e;

// 2. パラメータの宣言
parameters beta gamma varrho phi_pi phi_y rho_A rho_v;

// パラメータ値の代入
beta = 0.99;
gamma = 5;
varrho = 0.9;
phi_pi = 0.5;
phi_y = 0.5;
rho_A = 0.9;
rho_v = 0.7;

// 3. 方程式の定義
model(linear);
# kappa = (1-varrho)*(1-varrho*beta)*(gamma+1)/varrho;
x = x(+1)-(ii-ppi(+1))+(rho_A-1)*a;
ppi = beta*ppi(+1) + kappa*x + e;
ii = (1+phi_pi)*ppi + phi_y*x + v;
v = rho_v*v(-1) + z;
a = rho_A*a(-1) + epsilon;
end;

shocks;
var epsilon; stderr 0.5;
var z; stderr 0.5;
var e; stderr 0.5;
end;

// 全てゼロになることを念のためチェック
steady;

// モデルのチェック
check;

estimated_params;
gamma, gamma_pdf, 5, 2;
varrho, beta_pdf, 0.9, 0.05;
phi_pi, gamma_pdf, 0.5, 0.25;
phi_y, gamma_pdf, 0.5, 0.25;
rho_A, beta_pdf, 0.9, 0.05;
rho_v, beta_pdf, 0.7, 0.1;
stderr epsilon, inv_gamma_pdf, 5, 1;
stderr z, inv_gamma_pdf, 2, 1;
stderr e, inv_gamma_pdf, 1, 0.5;
end;

varobs x ppi ii;

estimation(datafile='dset.mat', mh_replic=125000, mh_drop = 0.2, mh_nblocks=2, mh_jscale=0.6, mode_compute = 4, mode_check);

chain1 = load(strcat('./NK_Linear_EST/metropolis/NK_Linear_EST_mh1_blck1.mat'))
chain2 = load(strcat('./NK_Linear_EST/metropolis/NK_Linear_EST_mh1_blck2.mat'));
csvwrite('chain1.csv',chain1.x2, 1, 0);
csvwrite('chain2.csv',chain2.x2, 1, 0);

