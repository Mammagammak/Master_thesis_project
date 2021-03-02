function difference_merged_neighbor_neurons_compare(diff_values)
% this function merges the interested *.mat files for further analysis
directory_to_save=input('enter the name of the folder you wish to save the current file to: e.g. "Conversion rates medium sensitivity NB w patch merged". = ');
oldFolder = cd(diff_values);
interval_spike_counts=input('Enter the spike counts of the interval varying experiment you want to analyze as a numeric variable: (if the experiment has standard protcol enter 0)=');
string=num2str(interval_spike_counts);
if interval_spike_counts~=0
cd(string);
else
string='_';
end

fList = dir('*.mat'); 

fList = {fList.name}';

for i=1:size(fList,1)
m=matfile(fList{i});


ultimate_max_pixel_values(i)=m.ultimate_max; %ultimate max pixel value in that experiment
ultimate_min_pixel_values(i)=m.ultimate_min; %ultimate min pixel value in that experiment

first_green_max_values(i)=m.first_green_max; %first green image max pixel value
first_green_min_values(i)=m.first_green_min; %first green image min pixel value

first_red_max_values(i)=m.first_red_max; %first red image max pixel value
first_red_min_values(i)=m.first_red_min; %first red image min pixel value

first_green_ROI_max_pixel_values(i)=m.first_green_ROI_max_pixel; %first green image max pixel value of the ROI
green_ROI_initial_value_mean_values(i)=m.green_ROI_initial_value_mean; %initial average brightness of the recorded neuron in that experiment

red_all_ROI_max_values(i)=m.red_all_ROI_max; % max red pixel value among all conditions within the ROI
red_all_ROI_min_values(i)=m.red_all_ROI_min; % min red pixel value among all conditions within the ROI

Rmax_values(i)=m.Rmax; %greatest red fluorescent condition's red fluorescent value
Rmax_index_values(i)=m.Rmax_index; % condition where red fluorescent was the greatest

normalization_pixel_values(i)=m.normalization_pixel; % pixel used to normalize that experiments values
normalized_green_initial_values(i)=(green_ROI_initial_value_mean_values(i)-ultimate_min_pixel_values(i))./normalization_pixel_values(i);

all_ROI_array{i}=m.ROI_array; % ROI from each experiment, this matrice contains 1's as the ROI and 0's as the other regions in the image
a{i}=find(all_ROI_array{i}{1}==1); %extract the ROI 
size_recorded(i)=size(a{i},1); %size of the recorded neuron as pixel count from each experiment

all_diff_total_values{i}=m.all_diff_values; %all conditions of all neurons from each experiment

Red_pre_diff_total_values{i}=all_diff_total_values{i}{3}; 
Red_pre_N_diff_total_values{i}=all_diff_total_values{i}{9};
Red_pre_N2_diff_total_values{i}=all_diff_total_values{i}{15};
Red_pre_N3_diff_total_values{i}=all_diff_total_values{i}{21}; % red pre condition conversion values from each experiment


clearvars m;
end    
    
% adding zeros in the beginning as the pre-pre comparison condition
 
for j=1:size(Red_pre_diff_total_values,2)
    
    carrier=Red_pre_diff_total_values{j};
    carrier=[0,carrier];
    Red_pre_diff_total_values_zero_added{j}=carrier;
    
    carrier_1=Red_pre_N_diff_total_values{j};
    carrier_1=[0,carrier_1]; 
    Red_pre_N_diff_total_values_zero_added{j}=carrier_1;
    
    carrier_2=Red_pre_N2_diff_total_values{j};
    carrier_2=[0,carrier_2];
    Red_pre_2_diff_total_values_zero_added{j}=carrier_2;
    
    carrier_3=Red_pre_N3_diff_total_values{j};
    carrier_3=[0,carrier_3];
    Red_pre_3_diff_total_values_zero_added{j}=carrier_3;    
end   

%finding the condition which the recorded neuron reaches the max
%conversion, then taking the same condition's corresponding conversion value at the non recorded neurons
for i=1:size(Red_pre_diff_total_values,2)
    conversion_max=max(Red_pre_diff_total_values{i});
    conversion_max_idx_logic=(Red_pre_diff_total_values{i}==conversion_max);
    conversion_max_idx=find(conversion_max_idx_logic==1);
    Conversion_final{1,i}(1)=Red_pre_diff_total_values{i}(conversion_max_idx);
    Conversion_final{1,i}(2)=Red_pre_N_diff_total_values{i}(conversion_max_idx);
    Conversion_final{1,i}(3)=Red_pre_N2_diff_total_values{i}(conversion_max_idx);
    Conversion_final{1,i}(4)=Red_pre_N3_diff_total_values{i}(conversion_max_idx);
end
Conversion_final_650=Conversion_final;
cd(oldFolder);
oldFolder=cd('Neighbor neurons compared high NB 100percent');
load('Conversion mergeddifference pixel wise no nucleus non recordeds included high sensitivity NB 100percent.mat')

f1=figure;
xlabel('Neurons');
ylabel('Change from the initial red fluorescence');
neurons=[1 2 3 4];
group_carrier=[1 2];

group=[];
y=[];
for i=1:size(Conversion_final,2)
    plot(neurons,Conversion_final{1,i},'color',rand(1,3),'LineStyle',':','Marker','.','MarkerSize', 30);
    y=[y,Conversion_final{i}(1,1),Conversion_final{i}(1,2)];
    group=[group,group_carrier];
   hold on
end


neurons=[1 2 3 4];
for i=1:size(Conversion_final_650,2)
    plot(neurons,Conversion_final_650{1,i},'color',rand(1,3),'LineStyle',':','Marker','.','MarkerSize', 30);
    group=[group,group_carrier];
    y=[y,Conversion_final_650{1,i}(1,1),Conversion_final_650{1,i}(1,2)];
   hold on
end
 hold off   

 [p,tbl,stats] = anova1(y,group,'on');


%% saving variables
cd(oldFolder);
cd(directory_to_save); 
s1='Conversion merged';
filename=strcat(s1,diff_values,string);
print(f1,filename,'-dsvg','-r600');

save(['Conversion 100 percent merged 650 and 850' diff_values '.mat']);
cd(oldFolder);

