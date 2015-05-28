function [data_bc, data_tc] = BAF_tQN(data_bc, data_tc)

flag1 = 0;
flag2 = 0;
if size(data_bc,2) > 1
    data_bc = data_bc';
    flag1 = 1;
end
if size(data_tc,2) > 1
    data_tc = data_tc';
    flag2 = 1;
end

data_bc_tQN = data_bc;
data_tc_tQN = data_tc;
data_baf = data_bc./data_tc;
tv_m = data_baf > 0 & data_baf < 1;
data_bc = data_bc(tv_m);
data_tc = data_tc(tv_m);

thre = 1.5;%0.9
data_ac = data_tc-data_bc;
a_fre = data_ac./data_tc;
b_fre = data_bc./data_tc;
[temp1,indx1] = sort(a_fre);
[temp2,indx2] = sort(b_fre);

temp = [temp1 temp2];
mean_values = mean(temp,2);
temp1 = zeros(length(a_fre),1);
temp2 = zeros(length(b_fre),1);
temp1(indx1) = mean_values;
temp2(indx2) = mean_values;
tv = temp1./a_fre > thre;
temp1(tv) = thre*a_fre(tv);
tv = temp2./b_fre > thre;
temp2(tv) = thre*b_fre(tv);
a_fre = temp1;
b_fre = temp2;

baf = b_fre./(a_fre+b_fre);
tv = baf > 0 & baf < 1;
d = round((data_tc(tv).*baf(tv)-data_bc(tv))./(1-baf(tv)));
data_bc(tv) = data_bc(tv)+d;
data_tc(tv) = data_tc(tv)+d;

data_bc_tQN(tv_m) = data_bc;
data_tc_tQN(tv_m) = data_tc;
data_bc = data_bc_tQN;
data_tc = data_tc_tQN;



if flag1
    data_bc = data_bc';
end
if flag2
    data_tc = data_tc';
end

end