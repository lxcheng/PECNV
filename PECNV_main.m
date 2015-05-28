function [PECNV_paras] = ...
    PECNV_main(init_PECNV_paras,depend_table,thres_EM,fn_nosuffix)
global mc_w
mc_w = 0.8;
thres_EM = 1e-4;
ds_info = [];
    thres_del = 0.006;
    init_PECNV_paras = PECNV_PECNV_paras(init_PECNV_paras);
    [PECNV_paras] = PECNV_screening...
        (init_PECNV_paras,depend_table,thres_EM,10);
end


function PECNV_paras = PECNV_PECNV_paras(init_PECNV_paras)
global tumor_range
PECNV_paras = cell(1,4);
%parameter initialization
    %---w---, ###variable in grid searching###
    if isempty(init_PECNV_paras{2}) %B w
           w_all = [0.1 0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9 0.8 0.9];
           o_all = [0.1 -0.2 -0.15 0 -0.2 0 0.1 -0.1 -0.05 0.1 0.05];
            % The soultions candidate are chosen by the tumor proportion range
            tv=w_all>=(1-tumor_range(2))&w_all<=(1-tumor_range(1));
            if sum(tv)==0
                %if no solution candidate, find the closest one, maybe two
                tv=(abs(w_all-(1-tumor_range(2)))==min(abs(w_all-(1-tumor_range(2)))));
                tv=tv|(abs(w_all-(1-tumor_range(2)))==min(abs(w_all-(1-tumor_range(1)))));
            end
            w_all=[w_all(tv) 1-(tumor_range(1)+tumor_range(2))/2];
            o_all=[o_all(tv) 0];
            clear tv
            w_0 = w_all;
    else %only one fixed (w,lrr) pair from previous searching
        w_0 = init_PECNV_paras{2};
    end
    N = size(w_0,2);
    PECNV_paras{2} = mat2cell(w_0,size(w_0,1),ones(1,size(w_0,2)));

    %---o---
    if isempty(init_PECNV_paras{1}) %B o
        o_0 = o_all;
    elseif length(init_PECNV_paras{1})==1%fixed lrr baseshift
        o_0 = repmat(init_PECNV_paras{1},1,N);
    else % multiple (w,lrr) pairs
        o_0 = init_PECNV_paras{1};
    end
    PECNV_paras{1} = mat2cell(o_0,size(o_0,1),ones(1,size(o_0,2)));
    %---varl---
    if isempty(init_PECNV_paras{3}) %B varl
        varl_0 = 0.2^2;
    else
        varl_0 = init_PECNV_paras{3};
    end
    PECNV_paras{3} = repmat({varl_0},1,N);

    %---varb---
    if isempty(init_PECNV_paras{4}) %B varb
        varb_0 = [(0.05^2)];
    else
        varb_0 = init_PECNV_paras{4};
    end
    PECNV_paras{4} = repmat({varb_0},1,N);

end
