function [PECNV_paras] = PECNV_screening(init_PECNV_paras,depend_table,thres1,max_iter1)

%---------------------run the algorithm------------------------------
%1xN cell vectors
o_all = init_PECNV_paras{1};
w_all = init_PECNV_paras{2};
varl_all = init_PECNV_paras{3};
varb_all = init_PECNV_paras{4};
PECNV_paras = cell(1,4);
t_all = 0; 
global gamma_sep
global gamma_sep_all;
gamma_sep_all=[];
clear init_PECNV_paras;init_PECNV_paras=cell(1,4);
for i=1:length(o_all)
    init_PECNV_paras(1) = o_all(i);
    init_PECNV_paras(3) = varl_all(i);
    init_PECNV_paras(4) = varb_all(i);
    init_PECNV_paras(2) = {w_all{i}(1)};
    gamma_sep=[];
    [o,w,varl,varb,nrIterations] = PECNV_EM_Newton_single_clone(init_PECNV_paras,depend_table,thres1,max_iter1);    
   gamma_sep_all=[gamma_sep_all {gamma_sep}];
    PECNV_paras{1} = [PECNV_paras{1} {o}];
    PECNV_paras{2} = [PECNV_paras{2} {w}];
    PECNV_paras{3} = [PECNV_paras{3} {varl}];
    PECNV_paras{4} = [PECNV_paras{4} {varb}];
    
        t = toc;
      %  disp('--------------- screening report -----------------')
        t_all = t_all+t;
end
end