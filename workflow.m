%excute the whole preprocess including preprocess,exomeCBS and mapprior;
%-----------------------------------------------------------------
%parameter that need to specify;
function workflow(namelist,out_path,GC_name,map_name)
%-------------------------------------------------------------------
fid_namelist=fopen(namelist);
if fid_namelist==-1
    error('can not find namelist!!!\n');
end
name_list=textscan(fid_namelist,'%s%s');
fclose(fid_namelist);
normalcount=name_list{1}{1};
normaldepth=name_list{2}{1};
for i=3:length(name_list{1})
    tumorcount=name_list{1}{i};
    tumordepth=name_list{2}{i};
    name=['s',num2str(i-1)];
    if exist([out_path,'\exome'],'dir')==0;
           mkdir([out_path,'\exome']);
    end
     if exist([out_path,'\seg'],'dir')==0;
           mkdir([out_path,'\seg']);
    end
outname=[out_path,'\',name,'.mat'];
[filterdata] = preprocess(normalcount,normaldepth,tumorcount,tumordepth,GC_name,map_name,outname);

segname=[out_path,'\seg',name,'.txt'];
matname=outname;
exomeCBS(segname,matname)

outexome=[out_path,'\exome','\',name,'exome.txt'];
outsegname=[out_path,'\seg','\',name,'seg.txt'];
matdata=outname;
mapprior(outexome,outsegname,matdata,segname);
end
CBSdatasource=[out_path,'\seg'];
exomesource=[out_path,'\exome'];
Outputsource=[out_path,'\results'];
if exist(Outputsource,'dir')==0
    mkdir(Outputsource);
end
PECNV(CBSdatasource,exomesource,Outputsource)
end
   


