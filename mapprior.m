
 function mapprior(outexome,outsegname,matdata,segname)
 %flag 用于判断baf信号是否为噪声，噪声为1，正常为0;
 %噪声判断标准为|baf-mean|>1.2*std;
 %newbaf用于计算对噪声进行处理得到一个新的baf值，newbaf=mean+(baf-mean)/3;
 %outexome:output filename that contain LRR BAF Flag information;
 %outsegname:output filename that contain segment and prior information
 %in each segment;
 %matdata:input filename that contain LRR and BAF produced by preprocess;
 %segname:input filename that contain segment information produced by
 %exomeCBS；
 %prior:input filename that contain proir information in each SNP position
 %produced by ASCAT results;
 %name='SA030';
 %outexome=['C:\Users\lxcheng\Desktop\cnv\program\preprocess\new\data\',name,'exome.txt'];
 fid_exome=fopen(outexome,'w');
 fprintf(fid_exome,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Chr','Startpos','Endpos','LRR','Flag','newlrr','BAF','Flag','NewBaf');
 %outsegname=['C:\Users\lxcheng\Desktop\cnv\program\preprocess\new\data\',name,'prior.txt'];
 fid=fopen(outsegname,'w');clear outsegname;
 %matdata=['C:\Users\lxcheng\Desktop\cnv\program\preprocess\new\data\',name,'.mat'];
 load(matdata);clear matdata;
 %segname=['C:\Users\lxcheng\Desktop\cnv\program\preprocess\new\data\','seg',name,'.txt'];
 fid_seg=fopen(segname);clear segname;
 temp_segdata=textscan(fid_seg,'%d%d%d%d','HeaderLines',1);
 segdata=[temp_segdata{1},temp_segdata{2},temp_segdata{3},temp_segdata{4}];
 clear temp_segdata;
% fid_prior=fopen('C:\Users\lxcheng\Desktop\cnv\program\preprocess\TNBCprior\priorprob.txt');
%  fid_prior=fopen(priorname);
%  temp_prior=textscan(fid_prior,'%f%f%f%f%f','HeaderLines',1);
%  prior=[temp_prior{1},temp_prior{2},temp_prior{3},temp_prior{4},temp_prior{5}];
%  clear temp_prior;
 %prior=priorprob;clear priorprob;
 %load('.\data\segdata.mat');
%  segdata_prior=zeros(length(segdata),6);
 %Chr startpos endpos prior_amlification prior_normal prior_deletion
%  segdata_prior(:,1:3)=segdata(:,1:3);
% %  templist=segdata_prior(:,6)<0.6;
% %  segdata_prior(templist,6)=0.5;
%  clear templist;
%  for i=1:length(segdata)
% %      if i==1 ||segdata(i,1)~=segdata(i-1,1)
% %            list=prior(:,1)==segdata(i,1);
% %            tempprior=prior(list,:);
% %            clear list;
% %      end
% %            list1=tempprior(:,2)>=segdata(i,2);
% %            list2=tempprior(:,2)<=segdata(i,3);
% %            list=list1&list2;
% %            newlist=find(list==1);
% %            if isempty(newlist)==1 || length(newlist)==1
% %                if sum(list1)==0
% %                    b=find(list2==1);
% %                    tempstart=b(length(b))-1;
% %                    tempend=b(length(b));
% %                else
% %                a=find(list1==1);
% %                tempstart=a(1);
% %                tempend=a(2);
% %                end
% %            else
% %            tempstart=newlist(1);
% %            tempend=newlist(length(newlist));
% %            end
% %            clear list list1 list2 newlist
%      
% %            tempstart=find(tempprior(:,2)>segdata(i,2));
% %            tempend=find(tempprior(:,2)>segdata(i,3));
% %           if length(tempstart)>1
% %               templist=prior(tempstart,1)==segdata(i,1);
% %               newtempstart=tempstart(templist);clear templist tempstart;
% %               tempstart=newtempstart;
% %               clear newtempstart;
% %           end
% %            if length(tempend)>1
% %               templist=prior(tempend,1)==segdata(i,1);
% %               newtempend=tempend(templist);clear templist tempend;
% %               tempend=newtempend;
% %               clear newtempend;
% %            end
%            area_a=0;area_n=0;area_d=0;
%            for j=1:(tempend-tempstart)
%                area_a=area_a+(tempprior(tempstart+j-1,3)+tempprior(tempstart+j,3))*(tempprior(tempstart+j,2)-tempprior(tempstart+j-1,2))*0.5;
%                area_n=area_n+(tempprior(tempstart+j-1,4)+tempprior(tempstart+j,4))*(tempprior(tempstart+j,2)-tempprior(tempstart+j-1,2))*0.5;
%                area_d=area_d+(tempprior(tempstart+j-1,5)+tempprior(tempstart+j,5))*(tempprior(tempstart+j,2)-tempprior(tempstart+j-1,2))*0.5;
%            end
%          segdata_prior(i,4)=area_a/(tempprior(tempend,2)-tempprior(tempstart,2));
%          segdata_prior(i,5)=area_n/(tempprior(tempend,2)-tempprior(tempstart,2));
%          segdata_prior(i,6)=area_d/(tempprior(tempend,2)-tempprior(tempstart,2));
%          segdata_prior(i,4)=segdata_prior(i,4)/(segdata_prior(i,4)+segdata_prior(i,5)+segdata_prior(i,6));
%          segdata_prior(i,5)=segdata_prior(i,5)/(segdata_prior(i,4)+segdata_prior(i,5)+segdata_prior(i,6));
%          segdata_prior(i,6)=segdata_prior(i,6)/(segdata_prior(i,4)+segdata_prior(i,5)+segdata_prior(i,6));
%  end
 %save('.\data\segdata_prior.mat','segdata_prior');
 fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n','Chr','startpos','endpos','prior_amplication','prior_normal','prior_deletion');
%  clear i;
 for i=1:length(segdata)-1
     fprintf(fid,'%d\t%d\t%d\t%f\t%f\t%f\n',segdata(i,1),segdata(i,2),segdata(i,3),0.33,0.33,0.33);
 end
 i=i+1;
 fprintf(fid,'%d\t%d\t%d\t%f\t%f\t%f',segdata(i,1),segdata(i,2),segdata(i,3),0.33,0.33,0.33);
 fclose(fid);
 
%-------------------------------------------------------------------------------------------------%
 chr=filterdata{1};
 pos=floor((double(filterdata{2})+double(filterdata{3}))/2);
 baf=filterdata{4};
 lrr=filterdata{5};
 flag=zeros(length(baf),1);
 newbaf=zeros(length(baf),1);
 lrrflag=zeros(length(lrr),1);
 newlrr=zeros(length(lrr),1);
 for i=1:length(segdata)
           list1=chr==segdata(i,1);
           list2=pos>=segdata(i,2);
           list3=pos<=segdata(i,3);
           list4=baf~=-1;
           list=list1&list2&list3&list4;
           tempbaf=baf(list);
           list_list=tempbaf<0.5;tempbaf(list_list)=1-tempbaf(list_list);
           clear list_list list4;
           tempflag=flag(list);
           tempnewbaf=newbaf(list);
          baf_mean=mean(tempbaf);
          baf_std=std(tempbaf,1);
          list4=abs(tempbaf-baf_mean)>1.2*baf_std;
          tempflag(list4)=1;
          tempnewbaf(list4)=(tempbaf(list4)-baf_mean)/3+baf_mean;
          a=find(list==1);
          flag(a)=tempflag;
          newbaf(a)=tempnewbaf;
          
          clear list list4 tempflag tempbaf tempnewbaf baf_mean baf_std;
          list=list1&list2&list3;
          templrr=lrr(list);
          templrrflag=lrrflag(list);
          tempnewlrr=newlrr(list);
          lrr_mean=mean(templrr);
          lrr_std=std(templrr,1);
          list4=abs(templrr-lrr_mean)>1.2*lrr_std;
          templrrflag(list4)=1;
          tempnewlrr(list4)=(templrr(list4)-lrr_mean)/3+lrr_mean;
          b=find(list==1);
          lrrflag(b)=templrrflag;
          newlrr(b)=tempnewlrr;
          clear a list1 list2 list3 list4 list b templrr templrrflag lrr_mean lrr_std ;
 end
 chrlist=unique(chr);
 for i=1:length(chrlist)-1
     list=chr==chrlist(i);
     tempchr=chr(list);
     tempstart=filterdata{2}(list);
     tempend=filterdata{3}(list);
     templrr=lrr(list);
     templrrflag=lrrflag(list);
     tempnewlrr=newlrr(list);
     tempbaf=baf(list);
     tempflag=flag(list);
     tempnewbaf=newbaf(list);
     for k=1:length(tempchr)
         fprintf(fid_exome,'%d\t%d\t%d\t%f\t%d\t%f\t%f\t%d\t%f\n',tempchr(k),tempstart(k),tempend(k),templrr(k),templrrflag(k),tempnewlrr(k),tempbaf(k),tempflag(k),tempnewbaf(k));
     end
     clear list tempchr tempstart tempend templrr tempbaf tempflag tempnewbaf;
 end
 i=i+1;
     list=chr==chrlist(i);
     tempchr=chr(list);
     tempstart=filterdata{2}(list);
     tempend=filterdata{3}(list);
     tempnewlrr=newlrr(list);
     templrr=lrr(list);
     templrrflag=lrrflag(list);
     tempbaf=baf(list);
     tempflag=flag(list);
     tempnewbaf=newbaf(list);
     for k=1:length(tempchr)-1
         fprintf(fid_exome,'%d\t%d\t%d\t%f\t%d\t%f\t%f\t%d\t%f\n',tempchr(k),tempstart(k),tempend(k),templrr(k),templrrflag(k),tempnewlrr(k),tempbaf(k),tempflag(k),tempnewbaf(k));
     end
     k=k+1;
     fprintf(fid_exome,'%d\t%d\t%d\t%f\t%d\t%f\t%f\t%d\t%f',tempchr(k),tempstart(k),tempend(k),templrr(k),templrrflag(k),tempnewlrr(k),tempbaf(k),tempflag(k),tempnewbaf(k));
 fclose(fid_exome);
 end
     
           
           
           