function threshold_sensitivity(diff_values)
% this function merges the interested *.mat files for further analysis
directory_to_save=input('enter the name of the folder you wish to save the current file to: e.g. "Conversion rates medium sensitivity NB w patch merged". = ');
oldFolder = cd(diff_values);
interval_spike_counts=input('Enter the spike counts of the interval you want to examine (if standard protcol enter 0)=');
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

% Polynom_values{i}=m.P; % coefficients of each experiment
%red_pre_polynom_coeffs{i}=Polynom_values{i}{3}; 
ultimate_max_pixel_values(i)=m.ultimate_max; %ultimate max pixel value in that experiment
ultimate_min_pixel_values(i)=m.ultimate_min; %ultimate max pixel value in that experiment

first_green_max_values(i)=m.first_green_max;
first_green_min_values(i)=m.first_green_min;

first_red_max_values(i)=m.first_red_max;
first_red_min_values(i)=m.first_red_min;

first_green_ROI_max_pixel_values(i)=m.first_green_ROI_max_pixel;
green_ROI_initial_value_mean_values(i)=m.green_ROI_initial_value_mean; %initial average brightness of the recorded neuron in that experiment

red_all_ROI_max_values(i)=m.red_all_ROI_max;
red_all_ROI_min_values(i)=m.red_all_ROI_min; 
Rmax_values(i)=m.Rmax; %greatest red fluorescent condition's red fluorescent value
Rmax_index_values(i)=m.Rmax_index; % condition where red fluorescent was the greatest

all_diff_total_values{i}=m.all_diff_values; %redpre condition of recorded neuron from each experiment
Red_pre_diff_total_values{i}=all_diff_total_values{i}{3};
normalization_pixel_values(i)=m.normalization_pixel; % pixel used to normalize that experiments values
normalized_green_initial_values(i)=(green_ROI_initial_value_mean_values(i)-ultimate_min_pixel_values(i))./normalization_pixel_values(i);

all_ROI_array{i}=m.ROI_array;
all_recorded_ROI_array{i}=all_ROI_array{i}{1};
a{i}=find(all_recorded_ROI_array{i}==1);
% all_ROI{i}=m.ROI;
% a{i}=find(all_ROI{i}==1);
b(i)=size(a{i},1);


Red_pre_diff_total_values_multiplied_size{i}=((Red_pre_diff_total_values{i}).*b(i));
Red_pre_diff_total_values_multiplied_size_divided_norm_greeninitial{i}=((Red_pre_diff_total_values{i}).*b(i)./((green_ROI_initial_value_mean_values(i)-ultimate_min_pixel_values(i))./normalization_pixel_values(i)));
Red_pre_diff_total_values_divided_norm_greeninitial{i}=(Red_pre_diff_total_values{i})./((green_ROI_initial_value_mean_values(i)-ultimate_min_pixel_values(i))./normalization_pixel_values(i));
Red_pre_diff_total_values_multiplied_size_divided_greeninitial{i}=((Red_pre_diff_total_values{i}).*b(i))./(green_ROI_initial_value_mean_values(i));
Red_pre_diff_total_values_divided_greeninitial{i}=(Red_pre_diff_total_values{i})./green_ROI_initial_value_mean_values(i);


%calculating global max of each exp and normalizing each experiment by its
%global maxima on multplied by size experiments
global_max_multiplied(i)=max(Red_pre_diff_total_values_multiplied_size{i});
global_min_multiplied(i)=min(Red_pre_diff_total_values_multiplied_size{i});
normalization_value_multiplied(i)=global_max_multiplied(i)-global_min_multiplied(i);
normalized_Red_pre_diff_total_values_multiplied_size{i}=((Red_pre_diff_total_values_multiplied_size{i}-global_min_multiplied(i)))./normalization_value_multiplied(i);

%calculating global max of each exp and normalizing each experiment by its
%global maxima 
global_max(i)=max(Red_pre_diff_total_values{i});
global_min(i)=min(Red_pre_diff_total_values{i});
normalization_value(i)=global_max(i)-global_min(i);
normalized_Red_pre_diff_total_values{i}=((Red_pre_diff_total_values{i}-global_min(i))./normalization_value(i));

summed_conversion_values(i)=sum(Red_pre_diff_total_values{i});
summed_norm_conversion_values(i)=sum(normalized_Red_pre_diff_total_values{i});

%linear regression for each experiment
Y=(normalized_Red_pre_diff_total_values{i})';
x=[10; 30; 80; 180; 430; 880; 1530; 2380];
X=[ones(size(x))  x];
[B,BINT,R,RINT,STATS]=regress(Y,X);
stats_regression(i)=STATS(1);
clearvars m;
 
end
% correlation analysis between conversion and the variables brightness and
% size
R_size_norm_conversion= corrcoef(summed_norm_conversion_values,b);
R_norm_brightness_norm_conversion=corrcoef(summed_norm_conversion_values,normalized_green_initial_values);
R_brightness_norm_conversion=corrcoef(summed_norm_conversion_values,green_ROI_initial_value_mean_values);

R_size_conversion=corrcoef(summed_conversion_values,b);
R_norm_brightness_conversion=corrcoef(summed_conversion_values,normalized_green_initial_values);
R_brightness_conversion=corrcoef(summed_conversion_values,green_ROI_initial_value_mean_values);

summed_conversion_values_tr=summed_conversion_values';
green_ROI_initial_value_mean_values_tr=green_ROI_initial_value_mean_values';
b_tr=b';

% X=[summed_conversion_values',green_ROI_initial_value_mean_values'];
% Y=[summed_conversion_values',b'];
% [R_brightness_conversion_plot,PValue_X] = corrplot(X);
% [R_size_conversion_plot,PValue_Y] = corrplot(Y);

% plot the regression line of the correlations
f8 = figure ('Color', [1 1 1]);
s1 = plot(green_ROI_initial_value_mean_values_tr, summed_conversion_values_tr, 'k+');
set(s1, 'MarkerSize', 8, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
set(l,'LineWidth', 2)
%%% axis display 
text(green_ROI_initial_value_mean_values_tr(end-1),summed_conversion_values_tr(end-1),num2str(R_brightness_conversion(1,2)),'Color','k');
xlabel('Green ROI initial brightness mean values', 'FontSize', 15)
ylabel('Summed conversion rates', 'FontSize', 15)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')

f9 = figure ('Color', [1 1 1]);
s2 = plot(b_tr,summed_conversion_values_tr, 'r+');
set(s2, 'MarkerSize', 8, 'LineWidth', 2);
%%% regression line
hold on
ll = lsline ;
set(ll,'LineWidth', 2)
%%% axis display 
text(b_tr(end-1),summed_conversion_values_tr(end-1),num2str(R_size_conversion(1,2)),'Color','k');
xlabel('Sizes of the neurons in pixels', 'FontSize', 15)
ylabel('Summed conversion rates', 'FontSize', 15)
set(gca, 'FontSize', 10, 'YMinorTick','on','XMinorTick','on')

% pca
merged_obs=[summed_conversion_values_tr,green_ROI_initial_value_mean_values_tr,b_tr];
[COEFF,SCORE,latent,tsquare]=princomp(zscore(merged_obs));


s=[];
if interval_spike_counts==0
    s=[10 30 80 180 430 880 1530 2380];
else 
    i=0;
    for k=1:5
    i=interval_spike_counts+i;
    s(k)=i;
    end
end
% for i=1:size(normalized_Red_pre_diff_total_values_multiplied_size,2)
%     threshold=(10:0.001:2380);
%     interp_y{i}=interp1(s,normalized_Red_pre_diff_total_values_multiplied_size{i},threshold, 'linear');
%     [I,J,V]=find(interp_y{i}==0.1);
%     threshold_spike(i)=(J+0.009);
%     
%     yData=normalized_Red_pre_diff_total_values_multiplied_size{i};
% plot( s,normalized_Red_pre_diff_total_values_multiplied_size{i},'m');
% hold on
% text(s(end-1),yData(end-1),num2str(b(i)),'Color','b');
% hold on
% text(s(end),yData(end),num2str(normalized_green_initial_values(i)),'Color','m');
% plot(threshold_spike,0.1, 'bp');
% end
% hold off
f1=figure;
figure(f1) % red pre diff raw values calculating by the max of the values who agrees the %10 and/or the %90 conditions
threshold_spikes_max=[];
saturation_spikes_max=[];
for i=1:size(normalized_Red_pre_diff_total_values_multiplied_size,2)  
    final_solution_max=0;
    final_solution_sat_max=0;
    for m=1:(size(s,2)-1)
    [P{i}{m}]=polyfit([s(m) s(m+1)],normalized_Red_pre_diff_total_values{i}(1,(m:m+1)),1);
    y_poly{m}=polyval(P{i}{m},[s(m) s(m+1)]);
    % threshold spikes==>solving the equation to find spike value corresponding to 10% of the max
    % conversion value
    % saturation spikes==>solving the equation to find spike value corresponding to 90% of the max
    % conversion value
    
    plot([s(m) s(m+1)],y_poly{m},'m');
    hold on
    
    syms x
    assume(x,'Real');
    assumeAlso(x>0);
    eq=P{i}{m}(1)*x+P{i}{m}(2)==0.1;
    eq2=P{i}{m}(1)*x+P{i}{m}(2)==0.9;
    solx_max=solve(eq,x);
    solx_sat_max=solve(eq2,x);
    solx_numerics_max=vpa(solx_max);
    solx_sat_numerics_max=vpa(solx_sat_max);
    candidate_solx_indice_max=[];
    candidate_solx_sat_indice_max=[];
    candidate_solx_indice_max=find(solx_numerics_max >= s(m) & solx_numerics_max <= s(m+1));
    candidate_solx_sat_indice_max=find(solx_sat_numerics_max >= s(m) & solx_sat_numerics_max <= s(m+1));
    
    if isempty(candidate_solx_indice_max) ~=0 
        candidate_solx_max=0;
    else
        for j=1:size(candidate_solx_indice_max,2)
        candidate_solx_max=solx_numerics_max(candidate_solx_indice_max(j));
        end
    end
    
      if isempty(candidate_solx_sat_indice_max) ~=0 
        candidate_solx_sat_max=0;
    else
        for j=1:size(candidate_solx_sat_indice_max,2)
        candidate_solx_sat_max=solx_sat_numerics_max(candidate_solx_sat_indice_max(j));
        end
    end
    max_solx_sat_max=max(candidate_solx_sat_max);
    max_solx_max=max(candidate_solx_max);
    if max_solx_max > final_solution_max
    final_solution_max=max_solx_max;
    end
    if max_solx_sat_max > final_solution_sat_max
    final_solution_sat_max=max_solx_sat_max;
    end
    end
%     text(s(7),y_poly(7),num2str(b(i)),'Color','b');
%     hold on
%     text(s(8),y_poly(8),num2str(normalized_green_initial_values(i)),'Color','k');
    threshold_spikes_max(i)=final_solution_max;
    saturation_spikes_max(i)=final_solution_sat_max;
    clearvars final_solution_max final_solution_sat_max;
end

% red pre diff raw values calculating by the min of the values who agrees the %10 and/or the %90 conditions
threshold_spikes_min=[];
saturation_spikes_min=[];
for i=1:size(normalized_Red_pre_diff_total_values_multiplied_size,2)  
    final_solution_min=0;
    final_solution_sat_min=0;
    for m=1:(size(s,2)-1)
    [P{i}{m}]=polyfit([s(m) s(m+1)],normalized_Red_pre_diff_total_values{i}(1,(m:m+1)),1);
    y_poly{m}=polyval(P{i}{m},[s(m) s(m+1)]);
    % threshold spikes==>solving the equation to find spike value corresponding to 10% of the max
    % conversion value
    % saturation spikes==>solving the equation to find spike value corresponding to 90% of the max
    % conversion value
    

    syms x
    assume(x,'Real');
    assumeAlso(x>0);
    eq=P{i}{m}(1)*x+P{i}{m}(2)==0.1;
    eq2=P{i}{m}(1)*x+P{i}{m}(2)==0.9;
    solx_min=solve(eq,x);
    solx_sat_min=solve(eq2,x);
    solx_numerics_min=vpa(solx_min);
    solx_sat_numerics_min=vpa(solx_sat_min);
    candidate_solx_indice_min=[];
    candidate_solx_sat_indice_min=[];
    candidate_solx_indice_min=find(solx_numerics_min >= s(m) & solx_numerics_min <= s(m+1));
    candidate_solx_sat_indice_min=find(solx_sat_numerics_min >= s(m) & solx_sat_numerics_min <= s(m+1));
    
    if isempty(candidate_solx_indice_min) ~=0 
        candidate_solx_min=0;
    else
        for j=1:size(candidate_solx_indice_min,2)
        candidate_solx_min=solx_numerics_min(candidate_solx_indice_min(j));
        end
    end
    
      if isempty(candidate_solx_sat_indice_min) ~=0 
        candidate_solx_sat_min=0;
    else
        for j=1:size(candidate_solx_sat_indice_min,2)
        candidate_solx_sat_min=solx_sat_numerics_min(candidate_solx_sat_indice_min(j));
        end
    end
    min_solx_sat_min=min(candidate_solx_sat_min);
    min_solx_min=min(candidate_solx_min);
    if min_solx_min > final_solution_min
    final_solution_min=min_solx_min;
    end
    if min_solx_sat_min > final_solution_sat_min
    final_solution_sat_min=min_solx_sat_min;
    end
    end
%     text(s(7),y_poly(7),num2str(b(i)),'Color','b');
%     hold on
%     text(s(8),y_poly(8),num2str(normalized_green_initial_values(i)),'Color','k');
    threshold_spikes_min(i)=final_solution_min;
    saturation_spikes_min(i)=final_solution_sat_min;
    clearvars final_solution_min final_solution_sat_min;
end
% calculating mean and standadr deviation 
threshold_spikes_min_mean=mean(threshold_spikes_min);
threshold_spikes_max_mean=mean(threshold_spikes_max);
saturation_spikes_min_mean=mean(saturation_spikes_min);
saturation_spikes_max_mean=mean(saturation_spikes_max);
threshold_spikes_min_std=std(threshold_spikes_min);
threshold_spikes_max_std=std(threshold_spikes_max);
saturation_spikes_min_std=std(saturation_spikes_min);
saturation_spikes_max_std=std(saturation_spikes_max);

%% saving variables
cd(oldFolder);
cd(directory_to_save); 
s1='Calculated with max of the threshold or saturation values';

filename1=strcat(s1,diff_values,string);


print(f1,filename1,'-dsvg','-r600');



save(['Threshold and saturation spike values' diff_values '.mat']);
cd(oldFolder);

