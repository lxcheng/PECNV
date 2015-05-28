fid=fopen('.\data\2_Tumor_file.txt');
data=textscan(fid,'%*s%d%d%f%f','HeaderLines',1,'TreatAsEmpty',{'NA','NAN'});
for i=1:10
    figure(i);
    templist=data{1}==i;
    scatter(data{2}(templist),data{3}(templist));
    clear templist;
end
fclose(fid);
clear;
