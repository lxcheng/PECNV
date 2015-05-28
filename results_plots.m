function results_plots(PECNV_paras,exomedata,depend_table,m,newdir)
%show results
% clear;close all;
   global gamma_sep_all;
   chr_indx=1;
   startpos_indx=2;
   endpos_indx=3;
   lrr_indx=4;
   flag_lrr=5;
   lrr_new=6;
   baf_indx=7;
   flag_baf=8;
   baf_new=9;
list_filter=(~exomedata{flag_lrr})&(~exomedata{flag_baf})&(exomedata{baf_indx}~=-1);
data_baf_sep = exomedata{baf_indx}(list_filter);
below_list=data_baf_sep<0.5;
data_baf_sep(below_list)=1-data_baf_sep(below_list);clear below_list;
data_lrr_sep = exomedata{lrr_indx}(list_filter);
data_chr_sep=exomedata{chr_indx}(list_filter);
data_pos_sep=floor((double(exomedata{startpos_indx}(list_filter))+double(exomedata{endpos_indx}(list_filter)))/2);
clear results temp exomedata;
w_all=PECNV_paras{1,2};
o_all=PECNV_paras{1,1};
for i=1:22
    close all;
    figure(i);
    [~,state_list]=max(gamma_sep_all{m});
    list=data_chr_sep==i;
    w=w_all{m};
    o=o_all{m};
    subplot(3,1,1);
    scatter(data_pos_sep(list),depend_table(state_list(list),3));
    ylabel('CN');
    subplot(3,1,2);hold on;
    scatter(data_pos_sep(list),data_baf_sep(list));
    y=w*2+depend_table(state_list(list),3)*(1-w);
    b=w+(1-w)*depend_table(state_list(list),3).*depend_table(state_list(list),4);
    scatter(data_pos_sep(list),b./y);
    ylabel('BAF');
    hold off;
    subplot(3,1,3);
    hold on;
    scatter(data_pos_sep(list),data_lrr_sep(list));
    yc=(1-w)*(depend_table(state_list(list),3))+w*2;
    scatter(data_pos_sep(list),log10(yc/2)+o);
    ylabel('Log Ratio');
    hold off;
    title(['chr',num2str(i)]);
    if exist([newdir,'\plots'],'dir')==0
        mkdir([newdir,'\plots']);
    end
    print(gcf,'-dpng',[newdir,'\plots\chr',num2str(i),'.png']);
    set(gcf,'visible','off');
end

