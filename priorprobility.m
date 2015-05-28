 clc;clear;  
  path='E:\files\graduation project\GPHMM_Release_Version_1.3\GPHMM_Release_Version_1.3\results\';
  fid3=fopen('.\data\2_Tumor_file.txt');
  info=textscan(fid3,'%*s%f%f%*f%*f','Headerlines',1,'TreatAsEmpty',{'na','nan','NA','NAN'});
  priorprob=zeros(length(info{1}),5);
  priorprob(:,1)=info{1}(:);priorprob(:,2)=info{2}(:);
  fclose(fid3);clear info;
  for i=1:112
      name=[path num2str(i) '_Tumor_file.results'];
      fid=fopen(name);
      data=textscan(fid,'%d%d%d%*d%*f%d%*d%*d%*s%*s%*d%*d','HeaderLines',13);
      fclose(fid);clear fid name;
      for j=1:length(data{1})
          tempstart=find(priorprob(:,2)==data{2}(j));
          tempend=find(priorprob(:,2)==data{3}(j));
          if length(tempstart)>1
              templist=priorprob(tempstart,1)==data{1}(j);
              newtempstart=tempstart(templist);clear templist tempstart;
              tempstart=newtempstart;
              clear newtempstart;
          end
           if length(tempend)>1
              templist=priorprob(tempend,1)==data{1}(j);
              newtempend=tempend(templist);clear templist tempend;
              tempend=newtempend;
              clear newtempend;
           end
          if data{4}(j)==2
              label=4;
          elseif data{4}(j)>2
              label=3;
          else
              label=5;
          end
          priorprob(tempstart:tempend,label)=priorprob(tempstart:tempend,label)+1;
      end
  end
  priorprob(:,3:5)=priorprob(:,3:5)./112;
  save('.\data\priorprob.mat','priorprob');
  outname='.\priorprob.txt';
  fid2=fopen(outname,'w');
  fprintf(fid2,'%s\t%s\t%s\t%s\t%s\n','Chr','Pos','amplification','normal','deletion');
  for i=1:length(priorprob)
      fprintf(fid2,'%d\t%d\t%d\t%d\t%d\n',priorprob(i,1),priorprob(i,2),priorprob(i,3),priorprob(i,4),priorprob(i,5));
  end
   fclose(fid2);
   clear;
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
          
      
      
      