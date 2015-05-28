function [bestindex,LL,acn]=findsolution(results,exomedata,PECNV_paras,state,depend_table,thres_del)
% name={'SA018','SA029','SA030','SA031','SA051'};
% num=4;thres_del=0.006;
% path=['C:\Users\lxcheng\Desktop\cnv\program\PECNV\newresults\',name{num},'\'];
% exomepath='C:\Users\lxcheng\Desktop\cnv\program\preprocess\pre_results\';
% priorname=[exomepath,name{num},'\',name{num},'prior.txt'];
% exomename=[exomepath,name{num},'\',name{num},'exome.txt'];
% fid=fopen(priorname);
% results = textscan(fid,'%d %d %d %f %f %f','HeaderLines',1);
% fclose(fid);
ampli=results{4};
normal=results{5};
deletion=results{6};
% load([path,'PECNV_paras']);
% load([path,'seginformation']);
% load([path,'state']);
%  fid_exome=fopen(exomename);
%    exomedata=textscan(fid_exome,'%d%d%d%f%d%f%f%d%f','HeaderLines',1);
%     fclose(fid_exome);
    chr_indx=1;
    startpos_indx=2;
    endpos_indx=3;
    lrr_indx=4;
    flag_lrr=5;
    lrr_new=6;
    baf_indx=7;
    flag_baf=8;        
    baf_new=9;
  
% depend_table = [...
%     1 1 0.01 0.5;...
%     2 1 1 1.0;...
%     3 1 2 0.5;...
%     4 1 2 1.0;...
%     5 1 3 0.67;...
%     6 1 3 1.0;...
%     7 1 4 0.75;...
%     8 1 4 0.5;...
%     9 1 4 1.0;...
%     10 1 5 0.8;...
%     11 1 5 0.6;...
%     12 1 5 1.0;...
%     13 1 6 5/6;...
%     14 1 6 4/6;...
%     15 1 6 0.5;...
%     16 1 6 1.0;...
%     17 1 7 6/7;...
%     18 1 7 5/7;...
%     19 1 7 4/7;...
%     20 1 7 1.0;...
%     ];
acn=zeros(length(state),1);
seginformation=results(1,1:3);
for i=1:length(state)
    statelist=state{i};
    a=sum(depend_table(statelist,3).*double(seginformation{3}-seginformation{2}+1));
    b=sum(seginformation{3}-seginformation{2}+1);
    acn(i)=a/double(b);
end
   list_filter=(~exomedata{flag_lrr})&(~exomedata{flag_baf})&(exomedata{baf_indx}~=-1);
   baf=exomedata{baf_indx}(list_filter);
  % baf(exomedata{flag_baf}==1)=exomedata{baf_new}(exomedata{flag_baf}==1);
   lrr=exomedata{lrr_indx}(list_filter);
  % lrr(exomedata{flag_lrr}==1)=exomedata{lrr_new}(exomedata{flag_lrr}==1);
   pos=floor((double(exomedata{startpos_indx}(list_filter))+double(exomedata{endpos_indx}(list_filter)))/2);  
   LL=zeros(length(state),1);total_del=zeros(length(state),1);
   for i=1:length(state)
       statelist=state{i};w=PECNV_paras{2}{i};o=PECNV_paras{1}{i};%varb=PECNV_paras{4}{i};
       varb=0;varl=0;num_lrr=0;num_baf=0;
       for k=1:length(seginformation{1})
           list_temp=exomedata{chr_indx}(list_filter)==seginformation{1}(k);
           lrr_temp=lrr(list_temp);
           baf_temp=baf(list_temp);
           pos_temp=pos(list_temp);
           list_temp_temp=pos_temp>=seginformation{2}(k)& pos_temp<=seginformation{3}(k);
           lrr_temp2=lrr_temp(list_temp_temp);
           baf_temp2=baf_temp(list_temp_temp);
           templist=baf_temp2>=0&baf_temp2<0.5;
           baf_temp2(templist)=1-baf_temp2(templist);
           w_all = w;%in single clone model, there is only one p
            ns = 2; %copy number of stromal cells
            mus = 0.5;% baf mean of stromal cells
           Y = w_all*ns+(1-w_all)*depend_table(statelist(k),3);
           Z = w_all*ns*mus+(1-w_all)*depend_table(statelist(k),3)*depend_table(statelist(k),4);
            LL(i)=sum((lrr_temp2-log10(Y/2)).^2)+LL(i);
            num_lrr=num_lrr+length(lrr_temp2);
            LL(i)=sum((baf_temp2(baf_temp2>=0)-Z/Y).^2)+LL(i);
            num_baf=num_baf+length(baf_temp2(baf_temp2>=0));
           if statelist(k)==0
               total_del(i)=total_del(i)+seginformation{3}-seginformation{2}+1;
           end
       end
       total_del(i)=double(total_del(i))/double(sum(seginformation{3}-seginformation{2}+1));
       LL(i)=LL(i)/num_lrr;
   end
   %--------------------------------------------------------------------%
   % find best index;
%   LL=LL.*num_lrr;
   list1=acn<4.51;
   list2=total_del<thres_del;
   list=list1&list2;
   a=1:12;
   b=a(list);
   if(sum(list)==0) bestindex=12;
   else
      [~,ii]=min(LL(list));
      bestindex=b(ii);
   end
end    
           
           
           
           
           
           
           
           
           
    
    
    
    
    
    
    
