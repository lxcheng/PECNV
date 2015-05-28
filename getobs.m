function [obslik] = getobs(w,o,varl,varb,data_baf,data_lrr,depend_table)
%function [obslik,condi_probs_fluct] = PECNV_get_obslik_single_clone(data_baf,data_lrr,o,w,varl,varb,depend_table)
%obslik:
%condi_probs
% magic_num = 0.02;
data_baf=data_baf';
data_lrr=data_lrr';
normal_prior=[1.5 1.5 1.5 1.5];
N = length(data_lrr); %number of data points
w_all = w;%in single clone model, there is only one p
ns = 2; %copy number of stromal cells
Nc = depend_table(:,3)'; %vector of copy numbers of different entries
mus = 0.5;% baf mean of stromal cells
Muc = depend_table(:,4)'; %vector of baf means of different entries
sigmal = sqrt(varl);
sigmab = sqrt(varb); 
Num_US = length(depend_table); % the number of unique states regulated by global parameters
tv_S = depend_table(:,2)~=0;
Y = w_all(depend_table(tv_S,2)')*ns+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S);
Z = w_all(depend_table(tv_S,2)')*ns*mus+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S).*Muc(tv_S);
US_indx = depend_table(tv_S,1);
obslik = zeros(sum(tv_S),N);
condi_probs_fluct = zeros(sum(tv_S&(depend_table(:,1)<=Num_US)),N);
for i=1:length(Y)
    if US_indx(i)<=Num_US %normal state
        %LRR
        obslik_LRR = PECNV_eval_pdf_LRR(data_lrr,log10(Y(i)/2)+o,sigmal);
        %BAF
        temp = find([1 3 8 15] == US_indx(i));
        if ~isempty(temp)
          obslik_BAF_Het = normal_prior(temp)*PECNV_eval_pdf_BAF(data_baf,Z(i)/Y(i),sigmab);
        %obslik_BAF = somatic_prob*obslik_BAF_Homo+(1-somatic_prob)*(Priors(1,:).*obslik_BAF_Homo+Priors(2,:).*obslik_BAF_Het);
        else
            obslik_BAF_Het = PECNV_eval_pdf_BAF(data_baf,Z(i)/Y(i),sigmab);
        end
        obslik_BAF=obslik_BAF_Het;
        if US_indx(i)==1
            fluct_prob = 0.01;
        else
            fluct_prob = 0.001;
        end
        if isempty(find(data_baf==-1))==1
           obslik(i,:) = (1-fluct_prob)*obslik_LRR.*obslik_BAF+fluct_prob*0.1;
        else 
           obslik(i,:) = (1-fluct_prob)*obslik_LRR+fluct_prob*0.1; 
        end
        condi_probs_fluct(i,:) = (fluct_prob*0.1)./obslik(i,:);
    else %the last state 
%         obslik(i,:) = PECNV_eval_pdf_LRR(data_lrr,-3.5,1.33).*PECNV_eval_pdf_BAF(data_baf,0.5,0.33);
%         obslik(i,:) = (1-fluct_prob)*PECNV_eval_pdf_LRR(data_lrr,-4,1).*PECNV_eval_pdf_BAF(data_baf,0.5,0.5/3)+fluct_prob*0.1;
    end
 end
end

%%%%%%%%%%%%%%%%%%%%%%%

function results = PECNV_eval_pdf_LRR(data,mu_lrr,Sigma_lrr)

if size(data,1)>size(data,2) %Nx1->1xN
    data = data';
end

if size(mu_lrr,1)>size(mu_lrr,2) %Nx1->1xN
    mu_lrr = mu_lrr';
end

if size(Sigma_lrr,1)>size(Sigma_lrr,2) %Nx1->1xN
    Sigma_lrr = Sigma_lrr';
end

results = normpdf(data,mu_lrr,Sigma_lrr);
end
%%%%%%%%%%%%%%%%%%%%%%%
function results = PECNV_eval_pdf_BAF(data,mu_baf,Sigma_baf)

if size(data,1)>size(data,2) %Nx1->1xN
    data = data';
end

if size(mu_baf,1)>size(mu_baf,2) %Nx1->1xN
    mu_baf = mu_baf';
end

if size(Sigma_baf,1)>size(Sigma_baf,2) %Nx1->1xN
    Sigma_baf = Sigma_baf';
end

results = normpdf(data,mu_baf,Sigma_baf);
end
