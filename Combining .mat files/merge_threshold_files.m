function merge_threshold_files(diff_values)
% this function merges the interested *.mat files for further analysis
directory_to_save=input('enter the name of the folder you wish to save the current file to: e.g. "Conversion rates medium sensitivity NB w patch merged". = ');
oldFolder = cd(diff_values);
n=input('Are you comparing different light intensity experiments?(if yes enter 1, if no enter 0)=');

fList = dir('*.mat'); 

fList = {fList.name}';

for i=1:size(fList,1)
m=matfile(fList{i});

threshold_spikes_min_all{i}=m.threshold_spikes_min;
saturation_spikes_min_all{i}=m.saturation_spikes_min;
threshold_spikes_max_all{i}=m.threshold_spikes_max;
saturation_spikes_max_all{i}=m.saturation_spikes_max;
end

% plotting distribution of threshold values based on sensitivity
f1=figure;
f2=figure;
f3=figure;
f4=figure;
figure(f1);
for i=1:size(fList,1)
h{i}=histfit(threshold_spikes_min_all{i},16,'kernel');
   if i==1& n==0
       a='High';
       b='w';
       c='k';
       d='k';
       e='--';
   elseif i==2 & n==0
       a='Low';
       b='w';
     c=[0.5 0.5 0.5];
       d=[0.5 0.5 0.5];
       e=':';
   elseif i==3 & n==0
       a='Medium';
       b='w';
       c='r';
       d='r';
       e='-';
   elseif i==1 & n==1
       a='High 100% light';
       b='w';
       c=[0 0 0.5];
       d=[0 0 0.5];
       e='-';
   elseif i==2 & n==1
       a='High 10% light';
       b='w';
       c='k';
       d='k';
       e='--';
   end
    h{i}(1).FaceColor=b;
    h{i}(1).FaceAlpha=0.3;
    h{i}(1).EdgeAlpha=0.4;
    h{i}(1).EdgeColor=d;
    h{i}(2).Color=c;
    h{i}(2).LineWidth=4;
    h{i}(2).LineStyle=e;
    h{i}(1).DisplayName=a;
    h{i}(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h{i}(2).DisplayName=a;
hold on
end
hold off

figure(f2);
for i=1:size(fList,1)
h{i}=histfit(threshold_spikes_max_all{i},16,'kernel');
   if i==1 & n==0
         a='High';
       b='w';
       c='k';
       d='k';
       e='--';
   elseif i==2 & n==0
       a='Low';
       b='w';
        c=[0.5 0.5 0.5];
       d=[0.5 0.5 0.5];
       e=':';
   elseif i==3 & n==0
       a='Medium';
       b='w';
       c='r';
       d='r';
       e='-';
    elseif i==1 & n==1
       a='High 100% light';
       b='w';
         c=[0 0 0.5];
       d=[0 0 0.5];
       e='-';
   elseif i==2 & n==1
       a='High 10% light';
       b='w';
       c='k';
       d='k';
       e='--';
   end
    h{i}(1).FaceColor=b;
    h{i}(1).FaceAlpha=0.3;
    h{i}(1).EdgeAlpha=0.4;
    h{i}(1).EdgeColor=d;
    h{i}(2).Color=c;
    h{i}(2).LineWidth=4;
    h{i}(2).LineStyle=e;
    h{i}(1).DisplayName=a;
    h{i}(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h{i}(2).DisplayName=a;
hold on
end
hold off

figure(f3);
for i=1:size(fList,1)
h{i}=histfit(saturation_spikes_min_all{i},16,'kernel');
    if i==1 & n==0
       a='High';
       b='w';
       c='k';
       d='k';
       e='--';
   elseif i==2 & n==0
       a='Low';
       b='w';
        c=[0.5 0.5 0.5];
       d=[0.5 0.5 0.5];
       e=':';
   elseif i==3 & n==0
       a='Medium';
       b='w';
       c='r';
       d='r';
       e='-';
 elseif i==1 & n==1
       a='High 100% light';
       b='w';
      c=[0 0 0.5];
       d=[0 0 0.5];
       e='-';
   elseif i==2 & n==1
       a='High 10% light';
       b='w';
       c='k';
       d='k';
       e='--';
   end
  h{i}(1).FaceColor=b;
    h{i}(1).FaceAlpha=0.3;
    h{i}(1).EdgeAlpha=0.4;
    h{i}(1).EdgeColor=d;
    h{i}(2).Color=c;
    h{i}(2).LineWidth=4;
    h{i}(2).LineStyle=e;
    h{i}(1).DisplayName=a;
    h{i}(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h{i}(2).DisplayName=a;
hold on
end
hold off

figure(f4);
for i=1:size(fList,1)
h{i}=histfit(saturation_spikes_max_all{i},16,'kernel');
   if i==1 & n==0
         a='High';
       b='w';
       c='k';
       d='k';
       e='--';
   elseif i==2 & n==0
       a='Low';
       b='w';
       c=[0.5 0.5 0.5];
       d=[0.5 0.5 0.5];
       e=':';
   elseif i==3 & n==0
       a='Medium';
       b='w';
       c='r';
       d='r';
       e='-';
  elseif i==1 & n==1
       a='High 100% light';
       b='w';
     c=[0 0 0.5];
       d=[0 0 0.5];
       e='-';
   elseif i==2 & n==1
       a='High 10% light';
       b='w';
       c='k';
       d='k';
       e='--';
   end
    h{i}(1).FaceColor=b;
    h{i}(1).FaceAlpha=0.3;
    h{i}(1).EdgeAlpha=0.4;
    h{i}(1).EdgeColor=d;
    h{i}(2).Color=c;
    h{i}(2).LineWidth=4;
    h{i}(2).LineStyle=e;
    h{i}(1).DisplayName=a;
    h{i}(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h{i}(2).DisplayName=a;
    
hold on
end
hold off


%% saving variables
cd(oldFolder);
cd(directory_to_save); 

filename1='Threshold spikes min values distributions plotted';
filename2='Threshold spikes max values distributions plotted';
filename3='Saturation spikes min values distributions plotted';
filename4='Saturation spikes max values distributions plotted';

print(f1,filename1,'-dsvg','-r600');
print(f2,filename2,'-dsvg','-r600');
print(f3,filename3,'-dsvg','-r600');
print(f4,filename4,'-dsvg','-r600');
save(['Threshold and saturation spike values distributions.mat']);
cd(oldFolder);



