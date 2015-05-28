%clc;clear;
function exomeCBS(segname,matname)
%segname:output filename that contains segment information;
%matname:intput filename that contains LRR and BAF data;
fid=fopen(segname,'w');
load(matname);
pos=floor((filterdata{2}+filterdata{3})/2);
lrr_cbs=[double(filterdata{1}),double(pos),double(filterdata{5})];
list=filterdata{4}~=-1;
temp=filterdata{4}(list);
list2=temp<0.5;
temp(list2)=1-temp(list2);
clear list2;
bafpos=floor((filterdata{2}(list)+filterdata{3}(list))/2);
baf_cbs=[double(filterdata{1}(list)),double(bafpos),double(temp)];
clear filterdata list bafpos temp;
lrr_segment=cghcbs(lrr_cbs,'showplot',false);
baf_segment=cghcbs(baf_cbs,'showplot',false);
results=[];m=1;
% smooth the segment;
for i=1:22
    lrr_startpos=lrr_segment.SegmentData(i).Start;
    lrr_endpos=lrr_segment.SegmentData(i).End;
    baf_startpos=baf_segment.SegmentData(i).Start;
    baf_startpos(1)=lrr_startpos(1);
    baf_endpos=baf_segment.SegmentData(i).End;
    baf_endpos(length(baf_endpos))=lrr_endpos(length(lrr_endpos));
    baf_mean=baf_segment.SegmentData(i).Mean;
    lrr_mean=lrr_segment.SegmentData(i).Mean;
    for j=1:length(lrr_startpos)
        startid=sum(lrr_startpos(j)>=baf_startpos);
        endid=sum(lrr_endpos(j)>baf_startpos);
        if startid==endid
            results(m,1)=i;
            results(m,2)=lrr_startpos(j);
            results(m,3)=lrr_endpos(j);
            results(m,4)=lrr_endpos(j)-lrr_startpos(j)+1;
            results(m,5)=lrr_mean(j);
            results(m,6)=baf_mean(startid);
            m=m+1;
        else
            for n=1:(endid-startid+1)
                if n==1
                   results(m,1)=i;
                   results(m,2)=lrr_startpos(j);
                   results(m,3)=baf_endpos(startid);
                   results(m,4)=baf_endpos(startid)-lrr_startpos(j)+1;
                   results(m,5)=lrr_mean(j);
                   results(m,6)=baf_mean(startid);
                   m=m+1;
                elseif n==endid-startid+1 && baf_startpos(startid+n-1)~=lrr_endpos(j)
                   results(m,1)=i;
                   results(m,2)=baf_startpos(startid+n-1);
                   results(m,3)=lrr_endpos(j);
                   results(m,4)=-baf_startpos(startid+n-1)+lrr_endpos(j)+1;
                   results(m,5)=lrr_mean(j);
                   results(m,6)=baf_mean(startid+n-2);
                   m=m+1;
                else
                   results(m,1)=i;
                   results(m,2)=baf_startpos(startid+n-1);
                   results(m,3)=baf_startpos(startid+n);
                   results(m,4)=baf_startpos(startid+n)-baf_startpos(startid+n-1)+1;
                   results(m,5)=lrr_mean(j);
                   results(m,6)=baf_mean(startid+n-2);
                   m=m+1;
                end
            end
        end
    end
end
 %save('.\data\segdata.mat','segdata');
 fprintf(fid,'%s\t%s\t%s\t%s\n','Chr','startpos','endpos','length');
 for i=1:length(results)-1
     if results(i,4)>1
     fprintf(fid,'%d\t%d\t%d\t%d\n',results(i,1),results(i,2),results(i,3),results(i,4));
     end
 end
 i=i+1;
 if results(i,4)>1
     fprintf(fid,'%d\t%d\t%d\t%d',results(i,1),results(i,2),results(i,3),results(i,4));
 end
 clear results;
 fclose(fid);
end
            
            
            
            
            
            
        
        
 