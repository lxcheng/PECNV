%获取baf lrr原始信号
function [filterdata] = get_baf_lrr(normalcount,normaldepth,tumorcount,tumordepth,GC_name,map_name,outname)
%本函数预处理外显子数据，包括过滤，去GC，计算RPKM，tumor,normal数据取比值后取log;
%   normalcount:normalcount文件名；
%   normaldepth:normaldepth文件名；
%   GC:GC文件名；
%  filterdata:处理完成数据，共5列：chr startpos endpos baf lrr
fid_normalcount=fopen(normalcount);
fid_normaldepth=fopen(normaldepth);
fid_tumordepth=fopen(tumordepth);
fid_tumorcount=fopen(tumorcount);
fid_GC=fopen(GC_name);
fid_map=fopen(map_name);
     if fid_normalcount==0
         error('normal count file can not been found');
     end
     if fid_normaldepth==0
         error('normal depth file can not been found');
     end
     if fid_tumorcount==0
         error('tumor count file can not been found');
     end
     if fid_tumordepth==0
         error('tumor depth file can not been found');
     end
     if fid_GC==0
         error('GC file can not been found');
     end
     count_normal=textscan(fid_normalcount,'%d%f%f%f','HeaderLines',1);
     depth_normal=textscan(fid_normaldepth,'%d%d%f%f','HeaderLines',1);
     count_tumor=textscan(fid_tumorcount,'%d%f%f%f','HeaderLines',1);
     depth_tumor=textscan(fid_tumordepth,'%d%d%f%f','HeaderLines',1);
     fclose(fid_normalcount);fclose(fid_normaldepth);fclose(fid_tumorcount);fclose(fid_tumordepth);
     
     GC=textscan(fid_GC,'%d%d%d%f','HeaderLines',1);
     map=textscan(fid_map,'%d%d%d%f','HeaderLines',1);
     %只选取常染色体
     chromelist=count_normal{1}>=1&count_normal{1}<=22;
     %1,计算RPKM
     totalcount_normal=sum(count_normal{4}(chromelist));
     RPKM_normal=count_normal{4}(chromelist)*(1.0e+9)./(totalcount_normal*(count_normal{3}(chromelist)-count_normal{2}(chromelist)));
     totalcount_tumor=sum(count_tumor{4}(chromelist));
     RPKM_tumor=count_tumor{4}(chromelist)*(1e+9)./(totalcount_tumor*(count_tumor{3}(chromelist)-count_tumor{2}(chromelist)));
     clear totalcount_normal totalcount_tumor
     %GC and mappbility;
     int_gc = floor(GC{4}(chromelist)*100);
     int_map = floor(map{4}(chromelist)*100);
      ncount_correct=RPKM_normal;
      tcount_correct=RPKM_tumor;
      m_all_normal = median(RPKM_normal);
      m_all_tumor=median(RPKM_tumor);
      
for i = min(int_gc):max(int_gc)
    for j = min(int_map):max(int_map)
        tv = int_gc == i & int_map == j;
        if sum(tv) > 0
            m_normal = median(RPKM_normal(tv));
            ncount_correct(tv) = round(RPKM_normal(tv)*m_all_normal/(m_normal+eps));
            m_tumor=median(RPKM_tumor(tv));
            tcount_correct(tv)= round(RPKM_tumor(tv)*m_all_tumor/(m_tumor+eps));
        end
    end
end
clear i j tv m_normal m_tumor int_gc int_map RPKM_tumor RPKM_normal
chr=count_normal{1}(chromelist);
startpos=count_normal{2}(chromelist);
endpos=count_normal{3}(chromelist);
clear chromelist
normallist=find(ncount_correct<100 & ncount_correct>5);
tumorlist=find(tcount_correct<100 & tcount_correct>5);
list=intersect(normallist,tumorlist);
newncount_correct=ncount_correct(list);
newtcount_correct=tcount_correct(list);
newchr=chr(list);
newstartpos=startpos(list);
newendpos=endpos(list);
clear chr startpos endpos
  clear normallist tumorlist list ncount_correct tcount_correct
  mean_normal=mean(newncount_correct);
  std_normal=std(newncount_correct,1);
  mean_tumor=mean(newtcount_correct);
  std_tumor=std(newtcount_correct,1);
  normallist=abs((newncount_correct-mean_normal)/std_normal)<2;
  tumorlist=abs((newtcount_correct-mean_tumor)/std_tumor)<2;
  list=normallist&tumorlist;
  
  chr=newchr(list);
  startpos=newstartpos(list);
  endpos=newendpos(list);
  ncount_correct=newncount_correct(list);
  tcount_correct=newtcount_correct(list);
  lrr=log10(tcount_correct./ncount_correct);
  
  clear list mean_normal std_normal mean_tumor std_tumor normallist tumorlist
  [depth_normal{3},depth_normal{4}]=BAF_tQN(depth_normal{3},depth_normal{4});
  baf_normal=depth_normal{3}./depth_normal{4};
  heterlist=baf_normal<0.75&baf_normal>0.25;
  temppos=depth_normal{2}(heterlist);
  tempchr=depth_normal{1}(heterlist);
  clear heterlist baf_normal
  [~,normalheterlist,tumorheterlist]=intersect(temppos,depth_tumor{2});
  templist=depth_tumor{1}(tumorheterlist)==tempchr(normalheterlist);
  heterlist=tumorheterlist(templist);
  clear templist;
 clear tumorheterlist normalheterlist
  heterchr=depth_tumor{1}(heterlist);
  heterpos=depth_tumor{2}(heterlist);
  [depth_tumor{3}(heterlist),depth_tumor{4}(heterlist)]=BAF_tQN(depth_tumor{3}(heterlist),depth_tumor{4}(heterlist));
  heterbaf=depth_tumor{3}(heterlist)./depth_tumor{4}(heterlist);
%  heterbaf(heterbaf<0.5)=1-heterbaf(heterbaf<0.5);
  clear heterlist temppos tempchr ncount_correct tcount_correct newncount_correct newtcount_correct
   chrlist=unique(chr);depthchr=chr;depthpos=startpos;depthsnp=zeros(length(depthpos),1)-1;
   chrlist=sort(chrlist);
   for k=1:length(chrlist)
       i=chrlist(k);
       temppos=heterpos(heterchr==i);
       tempbaf=heterbaf(heterchr==i);
       tempstart=startpos(chr==i);
       tempend=endpos(chr==i);
       newlist=find(chr==i);
       tempindex=newlist(1)-1;
       for j=1:length(temppos)
           startlist=tempstart<=temppos(j);
           endlist=tempend>=temppos(j);
           list=startlist&endlist;
           if sum(list)==1
               a=find(list==1);
               depthsnp(tempindex+a)=tempbaf(j);
           end
       end
   end
%    a=find(depthsnp==-1);
%    b=find(depthsnp~=-1);
%    if isempty(a)==0
%       for k=1:length(a)
%           chazhi=abs(a(k)-b);
%           [~,d]=min(chazhi);
%           depthsnp(a(k))=depthsnp(b(d));
%           clear chazhi
%       end
%    end
%   clear a b 
      filterdata=cell(1,5);
      filterdata{1}=chr;
      filterdata{2}=startpos;
      filterdata{3}=endpos;
      filterdata{4}=depthsnp;
      filterdata{5}=lrr; 
     save(outname,'filterdata');
end

