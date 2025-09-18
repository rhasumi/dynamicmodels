//=======================================
// Chapter 3, Ramsey model with tax 1
//   modified on 2020/01/02
//=======================================

// 1. 内生変数、外生変数の宣言
var C K;
varexo tauc tauk g;

// 2. パラメータの宣言
parameters alpha beta delta At;

// パラメータ値の代入
alpha = 0.3;
beta = 0.99;
delta = 0.25;
At = 1.0;

// 3. 方程式の定義
model;
(1+tauc(+1))*C(+1)/(1+tauc)/C-beta*((1-tauk(+1))*alpha*At*K^(alpha-1)-delta+1) = 0;
K-At*K(-1)^alpha-(1-delta)*K(-1)+C+g = 0;
end;

// 4. 初期値の定義
initval;
C = 1;
K = 1;
tauc = 0;
tauk = 0;
g = 0.1;
end;

steady;

// 5. 定常状態の計算と終端値の定義
endval;
C = 1;
K = 1;
tauc = 0.1;
tauk = 0;
g = 0.1;
end;

steady;

// 6. モデルのチェック
check;

// 7. シミュレーションの実行
shocks;
var tauc; periods 1:9; values 0;
end;

//simul(periods=31);
perfect_foresight_setup(periods=31);
perfect_foresight_solver;
for i = 1:size(M_.endo_names, 1)
  assignin('base', string(M_.endo_names(i)), oo_.endo_simul(i,:)');
end

// グラフ描写
figure(1)
subplot(2,1,1)
plot(0:30, C(1:31))
subplot(2,1,2)
plot(0:30, K(1:31))

