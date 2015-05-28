% filter segdata that segment length==1;
name='SA030';
inputname=['C:\Users\lxcheng\Desktop\cnv\program\preprocess\new\data\','seg',name,'.txt'];
outputname=['C:\Users\lxcheng\Desktop\cnv\program\preprocess\new\data\','newseg',name,'.txt'];
fid_in=fopen(inputname);
data=textscan(fid_in,'%d%d%d%d','HeaderLines',1);
fid_out=fopen(outputname,'w');
fprintf(fid_out,'%s\t%s\t%s\t%s\n','Chr','Startpos','Endpos','Length');
list=find(data{4}>1);
for i=1:length(list)
    k=list(i);
    fprintf(fid_out,'%d\t%d\t%d\t%d\n',data{1}(k),data{2}(k),data{3}(k),data{4}(k));
end
fclose(fid_out);
clear;
