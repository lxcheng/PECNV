function [state] = detection( results,exomedata,~,PECNV_paras,depend_table)
    global gamma_sep_all;
    chr_indx=1;
    startpos_indx=2;
    endpos_indx=3;
    lrr_indx=4;
    flag_lrr_indx=5;
    lrr_new_indx=6;
    baf_indx=7;
    flag_baf_indx=8;
    baf_new_indx=9;
    
    lrr_all=exomedata{lrr_indx};
    baf_all=exomedata{baf_indx};
    chr_all=exomedata{chr_indx};
    pos_all=floor((double(exomedata{startpos_indx})+double(exomedata{endpos_indx}))/2);
    lrrflag_all=exomedata{flag_lrr_indx};
    bafflag_all=exomedata{flag_baf_indx};
    newbaf_all=exomedata{baf_new_indx};
    newlrr_all=exomedata{lrr_new_indx};
    
%     ampli_indx=4;
%     normal_indx=5;
%     deletion_indx=6;
    
    chr_seg=results{chr_indx};
    startpos_seg=results{startpos_indx};
    endpos_seg=results{endpos_indx};
%     ampli_seg=results{ampli_indx};
%     normal_seg=results{normal_indx};
%     deletion_seg=results{deletion_indx};
    
    state=cell(1,length(gamma_sep_all));
    for i=1:length(gamma_sep_all)
        state{i}=zeros(length(chr_seg),1);
    end
    clear i;
    o_all=PECNV_paras{1};
    w_all=PECNV_paras{2};
    varl_all=PECNV_paras{3};
    varb_all=PECNV_paras{4};
    lrr_filter=lrrflag_all==1;
    lrr_all(lrr_filter)=newlrr_all(lrr_filter);
    baf_filter=bafflag_all==1;
    baf_all(baf_filter)=newbaf_all(baf_filter);
    for m=1:length(gamma_sep_all)
        o=o_all{m};
        w=w_all{m};
        varl=varl_all{m};
        varb=varb_all{m};
       for i=1:length(chr_seg)
           list=(chr_all==chr_seg(i))&(pos_all>=startpos_seg(i))&(pos_all<=endpos_seg(i));
           if sum(list)==0
              error('wrong segment\n');
           end
           lrr_temp=lrr_all(list);
           baf_temp=baf_all(list);
           pos_temp=pos_all(list);
           list_temp=baf_temp==-1;
           if sum(list_temp)>=1
               baf_temp_temp=baf_temp(~list_temp);
               pos_temp_temp=pos_temp(~list_temp);
               pos_temp_absent=pos_temp(list_temp);
               baf_temp_absent=zeros(length(pos_temp_absent),1);
               if isempty(baf_temp_temp)~=1
                  for k=1:length(pos_temp_absent)
                      [~,index]=min(abs(pos_temp_absent(k)-pos_temp_temp));
                      baf_temp_absent(k)=baf_temp_temp(index);
                      clear index;
                  end
                  baf_temp(list_temp)=baf_temp_absent;
                  clear list_temp;
               end
           end
          list_temp=baf_temp<0.5&baf_temp>=0;
          baf_temp(list_temp)=1-baf_temp(list_temp);clear list_temp;
           [obslik]=getobs(w,o,varl,varb,baf_temp,lrr_temp,depend_table);
           obs_sum=sum(log10(obslik),2);
           [~,index]=max(obs_sum);
           state{m}(i)=index;
       end
    end
end
    
    


