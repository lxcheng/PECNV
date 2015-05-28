function [varl,varb]=calculate_var
global data_lrr_sep;
global data_baf_sep;
[seg_info]=seglist;
seg_varl=[];
seg_varb=[];m=1;
for i=1:length(seg_info)
    if seg_info(i,1)==0 || seg_info(i,2)==0 || seg_info(i,2)-seg_info(i,1)<2
        continue;
    end
    startindex=seg_info(i,1);
    endindex=seg_info(i,2);
    lrr_temp=data_lrr_sep(startindex:endindex);
    seg_varl(m,1)=var(lrr_temp,1);
    seg_varl(m,2)=endindex-startindex+1;
    baf_temp=data_baf_sep(startindex:endindex);
    seg_varb(m,1)=var(baf_temp,1);
    seg_varb(m,2)=endindex-startindex+1;
    m=m+1;
end
 varb=sum(seg_varb(:,1).*seg_varb(:,2))/sum(seg_varb(:,2));
 varl=sum(seg_varl(:,1).*seg_varl(:,2))/sum(seg_varl(:,2));
end





