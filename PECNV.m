
function PECNV(CBSdatasource,exomesource,Outputsource)
%CBSdatasource :data that has been segmented by cghcbs;
%sourceflag: 1: directory otherwise: list file£»
% CBSdatasource='C:\Users\lxcheng\Desktop\cnv\program\PECNV\CBS';
% exomename='C:\Users\lxcheng\Desktop\cnv\program\PECNV\exomedata\SA030exome.txt';
% Outputsource='C:\Users\lxcheng\Desktop\cnv\program\PECNV\newresults';
global current_version
current_version = '1.0';
%     if nargin<2
%         error(['Insufficient input parameters, Please check again! ' ...
%             'More details in example.txt.'] );
%     end
depend_table = [...
    1 1 0.01 0.5;...
    2 1 1 1.0;...
    3 1 2 0.5;...
    4 1 2 1.0;...
    5 1 3 0.67;...
    6 1 3 1.0;...
    7 1 4 0.75;...
    8 1 4 0.5;...
    9 1 4 1.0;...
    10 1 5 0.8;...
    11 1 5 0.6;...
    12 1 5 1.0;...
    13 1 6 5/6;...
    14 1 6 4/6;...
    15 1 6 0.5;...
    16 1 6 1.0;...
    17 1 7 6/7;...
    18 1 7 5/7;...
    19 1 7 4/7;...
    20 1 7 1.0;...
    ];
thres_EM = 1e-4;
init_PECNV_paras = [{[]},{[]},{[]},{[]}]; 
global data_lrr_sep
global data_baf_sep
global data_startpos_sep
global data_endpos_sep
global data_ampli_sep
global data_normal_sep
global data_deletion_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global tumor_range
global gamma_sep_all
global data_pos_sep
global data_chr_sep
global data_chr
gamma_sep_all=[];
   disp(['PECNV (version ' current_version ') is loading...'])
   datafilelist = dir(CBSdatasource);
    if length(datafilelist)<3
        error('No CBS-array data files in the directory!');
    else 
       disp('-----PECNV will analyse ALL data files in directory-----')
       filename=cell(1,(length(datafilelist)-2));
       for i=3:length(datafilelist)
            filename{i-2}=datafilelist(i).name;
       end
       tumor_range_table=[0.1*ones((length(datafilelist)-2),1) ones((length(datafilelist)-2),1)];
       clear i
    end
        t_all = 0;
        if length(filename)==1
            disp(['-----PECNV batch annotation starts now, ONE SNP-array data file is found-----']);
        elseif length(filename)>1
            disp(['-----PECNV batch annotation starts now, ' num2str(length(filename)) ' SNP-array data files are found-----']);
        end
        for fileindx = 1:length(filename)
     %  for fileindx=3:3
      fprintf('deal with sample %d:%s \n',fileindx,filename{fileindx});
            data_lrr_sep = [];
            data_baf_sep = [];
            gamma_sep = [];
            condi_probs_sep = [];
            condi_probs_fluct_sep = [];
            %record time cost
            tic
            tumor_range=tumor_range_table(fileindx,:);
            results = regexp(filename{fileindx},'^(.+)\.+.+','tokens','once');
            if isempty(results)
                fn_nosuffix = filename{fileindx};
            else
                fn_nosuffix = results{1};
                if ~isempty(strfind(fn_nosuffix,'.'))
                    fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
                end
            end
            clear results;
            %--------open result files --------------
%             fid2 = fopen([Outputsource '/' fn_nosuffix '.txt'],'w');
%             if fid2 == -1
%                 error(['Can not open result file for ' filename{fileindx}]);
%             end
            
            %--------------load Segdata prior and exome data--------------------
            fid = fopen([CBSdatasource '/' filename{fileindx}],'r');
            if fid == -1
                error(['Can not open segmetation file: ' filename{fileindx}]);
            end
           temp=fgetl(fid);clear temp;
           results = textscan(fid,'%d %d %d %f %f %f', 'treatAsEmpty', {'NA', 'na'});% only autosomes are considered here
           fclose(fid);
             ampli_indx=4;
             normal_indx=5;
             deletion_indx=6;
             kk=strfind(fn_nosuffix,'seg');
             temp_name=fn_nosuffix(1:kk-1);
            exomename=[exomesource,'\',temp_name,'exome','.txt'];
            fid_exome=fopen(exomename);
            exomedata=textscan(fid_exome,'%d%d%d%f%d%f%f%d%f','HeaderLines',1);
            fclose(fid_exome);
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
        %    list_filter=(exomedata{baf_indx}~=-1);
            data_chr=results{chr_indx};
            data_startpos_sep = results{startpos_indx};
            data_endpos_sep=results{endpos_indx};
            data_ampli_sep=results{ampli_indx};
            data_normal_sep=results{normal_indx};
            data_deletion_sep=results{deletion_indx};
            
            data_baf_sep = exomedata{baf_indx}(list_filter);
            below_list=data_baf_sep<0.5;
            data_baf_sep(below_list)=1-data_baf_sep(below_list);clear below_list;
            data_lrr_sep = exomedata{lrr_indx}(list_filter);
            data_chr_sep=exomedata{chr_indx}(list_filter);
            data_pos_sep=floor((double(exomedata{startpos_indx}(list_filter))+double(exomedata{endpos_indx}(list_filter)))/2);
            clear temp
            
            fprintf(1,'Total %d exmome region are loaded from data file "%s".\n',length(data_lrr_sep),filename{fileindx});
            if size(data_startpos_sep,1)>size(data_startpos_sep,2) %make sure it's 1byN
                data_startpos_sep = data_startpos_sep';
            end
            if size(data_endpos_sep,1)>size(data_endpos_sep,2) %make sure     it's 1byN
                data_endpos_sep = data_endpos_sep';
            end
            if size(data_baf_sep,1)>size(data_baf_sep,2) %make sure it's 1byN
                data_baf_sep = data_baf_sep';
            end
            if size(data_lrr_sep,1)>size(data_lrr_sep,2) %make sure it's 1byN
                data_lrr_sep = data_lrr_sep';
            end
            if size(data_ampli_sep,1)>size(data_ampli_sep,2) %make sure it's 1byN
                data_ampli_sep = data_ampli_sep';
            end
             if size(data_normal_sep,1)>size(data_normal_sep,2) %make sure it's 1byN
                data_normal_sep = data_normal_sep';
             end
             if size(data_deletion_sep,1)>size(data_deletion_sep,2) %make sure it's 1byN
                data_deletion_sep = data_deletion_sep';
             end
            %------------------ call PECNV --------------------
            [PECNV_paras] = ...
                PECNV_main(init_PECNV_paras,depend_table,thres_EM,fn_nosuffix);
   %---------------------------------------------------%
   %detect CNV left behind
   state=[];
   state=detection(results,exomedata,list_filter,PECNV_paras,depend_table);
         newdir=[Outputsource,'\',temp_name]; 
    if exist(newdir,'dir')==0
        mkdir(newdir);
    end
 %   seginformation=[results(chr_indx),results(startpos_indx),results(endpos_indx),results()];
   thres_del=0.006;
    [bestindex,LL,acn]=findsolution(results,exomedata,PECNV_paras,state,depend_table,thres_del);
    results_plots(PECNV_paras,exomedata,depend_table,bestindex,newdir);
%     save([newdir,'\','PECNV_paras.mat'],'PECNV_paras');
%     save([newdir,'\','gamma_sep_all.mat'],'gamma_sep_all');
%     save([newdir,'\','state.mat'],'state');
%     save([newdir,'\','seginformation.mat'],'results');
%     save([newdir,'\','bestindex.mat'],'bestindex');
%     save([newdir,'\','LL.mat'],'LL');
%     save([newdir,'\','acn.mat'],'acn');
%---------------------------------------------------------------%
%Êä³ö½á¹û;
fid_new=fopen([newdir,'\','copy_number.txt'],'w');
fprintf(fid_new,'%s\t%s\t%s\t%s\t%s\n','Chr','StartPos','EndPos','CN_States','CN');
for i=1:length(results{1})
    fprintf(fid_new,'%d\t%d\t%d\t%d\n',results{1}(i),results{2}(i),results{3}(i),state{bestindex}(i),depend_table(state{bestindex}(i)),3);
end
end %if sourceflag  == 1
end