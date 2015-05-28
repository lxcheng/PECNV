function [o,w,varl,varb,nrIterations] = ...
    PECNV_EM_Newton_single_clone(init_PECNV_paras,depend_table, thresh, max_iter)
converged = 0;
num_iter = 1;
o = init_PECNV_paras{1};
w = init_PECNV_paras{2};
max_iter=100;
varl = init_PECNV_paras{3};
varb = init_PECNV_paras{4};
%---------------------------------------------------%
%caculate varl and varb;

%[varl,varb]=calculate_var;

%----------------------------------------------------%
while (num_iter <= max_iter) && ~converged
    % perform EM algorithm
    [obslik,LL,o_u,w_u,varl_u,varb_u] = ...
        PECNV_compute_ess_single_clone(o,w,varl,varb,depend_table);
        %o = o_u;
        o=o_u;
        w = w_u;
       varb = varb_u;
        varl=varl_u; 
       % fprintf('iter:%d time LL:%d o:%f\tw:%f\tvarb:%f\tvarl:%f\n',num_iter,LL,o,w,varb,varl);
%         if num_iter>1
%     converged = em_converged_m(obslik,previous_obslik,thresh);
%         end
%     previous_obslik = obslik; 
    num_iter =  num_iter + 1; 
end
nrIterations = num_iter - 1;

function [obslik,LL,o_u,w_u,varl_u,varb_u] = ...
    PECNV_compute_ess_single_clone(o,w,varl,varb,depend_table)
global data_lrr_sep
global data_baf_sep
global gamma_sep
global condi_probs_fluct_sep
global tumor_range
global data_ampli_sep
global data_normal_sep
global data_deletion_sep
%-----------------------E step-----------------------------
N = 0; % the size of the whole data set
obs_baf = data_baf_sep;
obs_lrr = data_lrr_sep;
N = N+length(obs_lrr);
    %condi_probs: Pi(G=k|S=j,O)
    [obslik,condi_probs_fluct] = PECNV_get_obslik_single_clone(obs_baf,obs_lrr,o,w,varl,varb,depend_table);
    %     [alpha, beta, gamma, current_ll, xi_summed] = fwdback(prior, transmat, obslik);
   % [alpha, gamma, current_ll, beta, xi_summed] = Forward_Backward_Algorithm(prior, transmat, obslik);
    %         [path,current_ll] = viterbi_path(prior, transmat, obslik);
   % prior_all=[data_deletion_sep;data_deletion_sep;data_normal_sep;data_normal_sep;repmat(data_ampli_sep',1,16)'];
   % gamma_sep=obslik.*prior_all;
    [seg_info]=seglist;gamma_sep=obslik;
    for m=1:length(seg_info)
        index1=seg_info(m,1);
        index2=seg_info(m,2);
        if index1==0||index2==0
            continue;
        end
        temp_obslik=obslik(:,index1:index2);
        segnum=floor((index2-index1+1)/20);
        temp_prior=[data_deletion_sep(m);data_deletion_sep(m);data_normal_sep(m);data_normal_sep(m);repmat(data_ampli_sep(m)',1,length(depend_table)-4)'];
        if segnum>=1
            for t=1:segnum
                newtemp_obslik=temp_obslik(:,20*(t-1)+1:20*t);
                a=prod(newtemp_obslik,2);
                b=a/sum(a);
                c=b.*temp_prior;
                gamma_sep(:,index1+20*(t-1):index1+20*t-1)=repmat(c,1,20);
            end
        end
        if 20*segnum==index2-index1+1
            continue;
        else
        t=segnum+1;
                newtemp_obslik=temp_obslik(:,20*(t-1)+1:(index2-index1+1));
                a=prod(newtemp_obslik,2);
                b=a/sum(a);
                c=b.*temp_prior;
                gamma_sep(:,index1+20*(t-1):index2)=repmat(c,1,index2-index1-20*(t-1)+1);
        end
%         sum_sum=sum(temp_obslik,1);
%         for t=1:length(sum_sum)
%             temp_obslik(:,t)=temp_obslik(:,t)/sum_sum(t);
%         end
     %   clear sum_sum;
%         log_tempobslik=log10(temp_obslik);
%         sum_log=sum(log_tempobslik,2);
%         [~,templist]=sort(sum_log,'descend');rate_list=ones(length(templist)-1,1);
%         for c=1:length(templist)-1
%             small=temp_obslik(templist(c+1),:);
%             big=temp_obslik(templist(c),:);
%             rate_list(c)=prod(small./big,2);
%         end
%         clear c;
%         real_rate=ones(length(templist),1);
%         for c=2:length(real_rate)
%             real_rate(c)=real_rate(c-1)*rate_list(c-1);
%         end
%         clear c rate_list;
%         real_rate=real_rate/sum(real_rate);
%       %  real_rate(real_rate<1e-24)=1e-24;
%         mean_obs=ones(length(templist),1);
%         mean_obs(templist)=real_rate;clear real_rate templist;
%         mean_obs=mean(temp_obslik,2);
%         mean_obs=mean_obs/sum(mean_obs);
%        temp_prior=[data_deletion_sep(m);data_deletion_sep(m);data_normal_sep(m);data_normal_sep(m);repmat(data_ampli_sep(m)',1,length(depend_table)-4)'];
%       gamma_sep(:,index1:index2)=repmat(mean_obs.*temp_prior,1,index2-index1+1);
%         if isempty(isnan(gamma_sep(1,:)))==0
%              fprintf('m:%d\n',m);
%         end
        clear index1 index2 ;
    end
    %LL=sum(sum(gamma_sep));
    for k=1:length(gamma_sep)
        gamma_sep(:,k)=gamma_sep(:,k)./sum(gamma_sep(:,k));
    end
    [~,statelist]=max(gamma_sep);LL=0;
    for mm=1:length(statelist)
       LL=LL+log10(obslik(statelist(mm),mm));
    end
    clear mm;
    for tt=1:length(seg_info)
        if seg_info(tt,1)==0 || seg_info(tt,2)==0
            continue;
        end
        pp=floor(double(seg_info(tt,1)+seg_info(tt,2))/2);
        if statelist(pp)>=5
            LL=LL+log10(data_ampli_sep(tt));
        elseif statelist(pp)<=2
            LL=LL+log10(data_deletion_sep(tt));
        else 
            LL=LL+log10(data_normal_sep(tt));
        end
    end
    clear k;
    condi_probs_fluct_sep = condi_probs_fluct;
    clear condi_probs_fluct;
    

%-----------------------M step-----------------------------
%update global parameters
% update w
w_tol = 0.01;
max_iter = 30;
w_u = PECNV_update_w_single_clone(o,w,varl,varb,depend_table,w_tol,max_iter);
if w_u>(1-tumor_range(1)),w_u =(1-tumor_range(1)); end
if w_u<(1-tumor_range(2)),w_u = (1-tumor_range(2)); end

%update o
o_u = PECNV_update_o_single_clone(w_u,o,depend_table);
%varl_u=varl;
%varb_u=varb;

%update varl
varl_u = PECNV_update_varl_single_clone(o_u,w_u,depend_table);

%update varb
varb_u = PECNV_update_varb_single_clone(w_u,depend_table);


%--------------------------------------------------------------------------
function o_u = PECNV_update_o_single_clone(w,o,depend_table)
global data_lrr_sep
%global data_baf_sep
global gamma_sep
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
% mus = 0.5;% baf mean of stromal cells
% Muc = depend_table(:,4)'; %row vector of baf means of different entries
Num_US = length(depend_table); % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);

iter = 0;
while 1
    Y = w*ns+(1-w)*Nc(tv);
%     Z = w*ns*mus+(1-w)*Nc(tv).*Muc(tv);
%     Y_d = ns-Nc(tv);
%     Z_d = ns*mus-Nc(tv).*Muc(tv);
   
    %first order differential
    ELL_L_D_1 = 0;
    %second order differential
    ELL_L_D_2 = 0;
    %obs_baf = data_baf_sep;
    obs_lrr = data_lrr_sep;
    post_probs_not_fluct = gamma_sep(1:length(Y),:);
  %  post_probs = gamma_sep(1:length(Y),:);
    part1 = o-obs_lrr;        
        
        for i=1:length(Y) 
            %LRR
            %first order differential
            ELL_L_D_1 = ELL_L_D_1 + sum(post_probs_not_fluct(i,:).*((log10(Y(i)/2)+part1)));
            %second order differential
            ELL_L_D_2 = ELL_L_D_2 + sum(post_probs_not_fluct(i,:));            
        end
    

%     ELL_L_D_1 = -1/(log(10)*varl)*ELL_L_D_1;
%     ELL_B_D_1 = -1/(varb)*ELL_B_D_1;
%     ELL_ALL_D_1 = ELL_L_D_1+ELL_B_D_1;
% 
%     ELL_L_D_2 = -1/(log(10)*varl)*ELL_L_D_2;
%     ELL_B_D_2 = -1/(varb)*ELL_B_D_2;
%     ELL_ALL_D_2 = ELL_L_D_2+ELL_B_D_2;
   w_tol=0.01;max_iter=30;
    o_u = o- ELL_L_D_1/ELL_L_D_2;
    iter = iter+1;  
    if abs(o_u-o)<w_tol || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        o_u =o;
    end
end

%--------------------------------------------------------------------------
% function o_u = PECNV_update_o_single_clone(w,depend_table)
% global data_lrr_sep
% global gamma_sep
% ns = 2; %copy number of stromal cells
% Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
% Num_US = length(depend_table); % the number of unique states regulated by global parameters
% tv = (depend_table(:,1)<=Num_US);
% Y = w*ns+(1-w)*Nc(tv);
% 
% o_num = 0;
% o_den = 0;
% obs_lrr = data_lrr_sep;
%     post_probs_not_fluct = gamma_sep(1:length(Y),:);
%     for i=1:length(Y)
%         o_num = o_num + sum(post_probs_not_fluct(i,:).*(obs_lrr-(log10(Y(i)/2)*1)));
%     end
%     o_den = o_den + sum(sum(post_probs_not_fluct));
% o_u = o_num/o_den;
%--------------------------------------------------------------------------
function w_u = PECNV_update_w_single_clone(o,w,varl,varb,depend_table,w_tol,max_iter)
global data_lrr_sep
global gamma_sep
global data_baf_sep
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
mus = 0.5;% baf mean of stromal cells
Muc = depend_table(:,4)'; %row vector of baf means of different entries
Num_US = length(depend_table); % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);

iter = 0;
while 1
    Y = w*ns+(1-w)*Nc(tv);
    Z = w*ns*mus+(1-w)*Nc(tv).*Muc(tv);
    Y_d = ns-Nc(tv);
    Z_d = ns*mus-Nc(tv).*Muc(tv);
   
    %first order differential
    ELL_L_D_1 = 0;
    ELL_B_D_1 = 0;
    %second order differential
    ELL_L_D_2 = 0;
    ELL_B_D_2 = 0;
    
    obs_baf = data_baf_sep;
    obs_lrr = data_lrr_sep;
    post_probs_not_fluct = gamma_sep(1:length(Y),:);
    post_probs = gamma_sep(1:length(Y),:);
    part1 = o-obs_lrr;        
        
        for i=1:length(Y) 
            %LRR
            %first order differential
            ELL_L_D_1 = ELL_L_D_1 + sum(post_probs_not_fluct(i,:).*((log10(Y(i)/2)+part1)*Y_d(i)/Y(i)));
            %second order differential
            ELL_L_D_2 = ELL_L_D_2 + sum(post_probs_not_fluct(i,:).*((1/log(10)-(1*log10(Y(i)/2)+part1))*(Y_d(i)/Y(i))^2));
            %BAF
            %first order fifferential
            ELL_B_D_1 = ELL_B_D_1 + sum(post_probs(i,:).*((Z_d(i)*Y(i)-Z(i)*Y_d(i))*(Z(i)-obs_baf.*Y(i))/Y(i)^3));
            %second order differential
            ELL_B_D_2 = ELL_B_D_2 + sum(post_probs(i,:).*((Z_d(i)*Y(i)-Z(i)*Y_d(i))*(Z_d(i)*Y(i)+Y_d(i)*(2*Y(i)*obs_baf-3*Z(i)))/Y(i)^4));
        end
    

    ELL_L_D_1 = -1/(log(10)*varl)*ELL_L_D_1;
    ELL_B_D_1 = -1/(varb)*ELL_B_D_1;
    ELL_ALL_D_1 = ELL_L_D_1+ELL_B_D_1;

    ELL_L_D_2 = -1/(log(10)*varl)*ELL_L_D_2;
    ELL_B_D_2 = -1/(varb)*ELL_B_D_2;
    ELL_ALL_D_2 = ELL_L_D_2+ELL_B_D_2;

    w_u = w - ELL_ALL_D_1/ELL_ALL_D_2;
    iter = iter+1;
    
    if abs(w_u-w)<w_tol || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        w = w_u;
    end
end
%--------------------------------------------------------------------------
function varl_u = PECNV_update_varl_single_clone(o,w,depend_table)
global data_lrr_sep
global gamma_sep
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
Num_US = length(depend_table); % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Y = w*ns+(1-w)*Nc(tv);
varl_num = 0; %numerator of o estimator
varl_den = 0;
 obs_lrr = data_lrr_sep;
 post_probs_not_fluct = gamma_sep(1:length(Y),:);
    for i=1:length(Y)
        varl_num = varl_num +  sum(post_probs_not_fluct(i,:).*(obs_lrr-(log10(Y(i)/2)+o)).^2);
    end
    varl_den = varl_den + sum(sum(post_probs_not_fluct));
varl_u = varl_num/varl_den;

%--------------------------------------------------------------------------
function [varb_u] = PECNV_update_varb_single_clone(w,depend_table)
global data_baf_sep
global gamma_sep
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %row vector of copy numbers of different entries
mus = 0.5;% baf mean of stromal cells
Muc = depend_table(:,4)'; %row vector of baf means of different entries
Num_US = length(depend_table); % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Y = w*ns+(1-w)*Nc(tv);
Z = w*ns*mus+(1-w)*Nc(tv).*Muc(tv);
varb_num = 0;
varb_den = 0;
obs_baf = data_baf_sep;
    post_probs =  gamma_sep(1:length(Y),:);
    varb_den = varb_den + sum(sum(post_probs));
     for i=1:length(Y)
        varb_num = varb_num + sum(post_probs(i,:).*(obs_baf-Z(i)/Y(i)).^2);
     end         
varb_u = varb_num/varb_den;



