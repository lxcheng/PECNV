function [seg_info]=seglist
% get segmentation information list
  global data_chr
  global data_chr_sep
  global data_pos_sep
  global data_startpos_sep
  global data_endpos_sep
% fid_seg=fopen('C:\Users\lxcheng\Desktop\cnv\program\preprocess\segprior\SA018prior.txt');
% seg=textscan(fid_seg,'%d%d%d%*f%*f%*f','HeaderLines',1);
% data_chr=seg{1};
% data_startpos_sep=seg{2};
% data_endpos_sep=seg{3};
% clear seg;fclose(fid_seg);
% fid_exome=fopen('C:\Users\lxcheng\Desktop\cnv\program\preprocess\data\SA018exome.txt');
% exomedata=textscan(fid_exome,'%d%d%d%*f%*f','HeaderLines',1);
% data_chr_sep=exomedata{1};
% data_pos_sep=floor((double(exomedata{2})+double(exomedata{3}))/2);
% fclose(fid_exome);clear exomedata;
  seg_info=zeros(length(data_chr),2);
   for i=1:length(data_chr)
       chr_list=data_chr_sep==data_chr(i);
       startpos_list=data_pos_sep>=data_startpos_sep(i);
       endpos_list=data_pos_sep<=data_endpos_sep(i);
%        startindex=chr_list&startpos_list;
%        endindex=chr_list&endpos_list;
%        if sum(startindex)==0 || sum(endindex)==0
%            fprintf('wrong\n');
%            return;
%        end
%        newstartindex=find(startindex==1);
%        newendindex=find(endindex==1);
  list=chr_list&startpos_list&endpos_list;
  list_list=find(list==1);
  if sum(list_list)==0
      continue;
  end
       seg_info(i,:)=[list_list(1),list_list(length(list_list))];
   end
end

