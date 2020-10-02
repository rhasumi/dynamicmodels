%=======================================
% Chapter 8, Estimation
%   modified on 2020/01/02
%=======================================

% パスのセット
cd 'C:\Users\USERNAME\ch8_dynare_test'

% データセットの作成
script_dataset

% Dynareのパスのセット
addpath C:\programs\dynare\4.5.7\matlab

% パラメータ推定、出力
dynare NK_Linear_EST

