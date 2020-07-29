%=======================================
% Chapter 8, Preparing data
%   modified on 2019/10/19
%=======================================

FqReport8 = csvread('FqReport8.csv', 1, 1);

GAP0 = FqReport8(:,1);
PGDP0 = FqReport8(:,2);
IRS0 = FqReport8(:,3);

PC_PGDP0 = [NaN; (PGDP0(2:end)./PGDP0(1:end-1)-1)*100];

tt = 2:77;
MGAP = mean(GAP0(tt));
MPC_PGDP = mean(PC_PGDP0(tt));
MIRS = mean(IRS0(tt));

GAP = GAP0(tt)-MGAP;
PC_PGDP = PC_PGDP0(tt)-MPC_PGDP;
IRS4 = (IRS0(tt)-MIRS)*0.25;

x = GAP;
ppi = PC_PGDP;
ii = IRS4;

save dset.mat x ppi ii;

out2 = [GAP0(tt) PC_PGDP0(tt) IRS0(tt)];
lab2 = {'GAP', 'PC_PGDP', 'IRS'};
util_csvwrite('dset_fig.csv', out2, lab2);
