function difference_pixel_wise_no_nucleus_non_recorded_inc(before,after,tags,finalfile)
% calculates the conversion rate as the increase in percentage by comparing
% each condition with the pre condition (condition before any electrical stimulation)

%!! in this function, the conversion rate of the non recorded (non-patched) neurons 
%are also calculated.

% INPUTS: 
% BEFORE and AFTER are structural arrays with at least one string, name of 
% the image files.  
% TAGS is a structural array with legends to be used for the images.
% FINALFILE is the name you would like to catalog the data under.  This is
% typically the name of the experiment. 
%
% Sample entry
% difference_pixel_wise_no_nucleus_non_recorded_inc({'160215_E7_Preconversion_Green.tif';'160215_E7_Preconversion_Red.tif'},{'160215_E7_Postconversion_Green_Pulsed light.tif';'160215_E7_Postconversion_Red_Pulsed light.tif'},{'Green';'Red'},'160215')


if (nargin<1)
before = {'180608_E1_img_000000001_Default_000.tif';'180608_E1_img_000000002_Default_000.tif'};
after={'180608_E1_img_000000004_Default_000.tif';'180608_E1_img_000000003_Default_000.tif';'180608_E1_img_000000006_Default_000.tif';'180608_E1_img_000000005_Default_000.tif';'180608_E1_img_000000008_Default_000.tif';'180608_E1_img_000000007_Default_000.tif';'180608_E1_img_000000010_Default_000.tif';'180608_E1_img_000000009_Default_000.tif';'180608_E1_img_000000012_Default_000.tif';'180608_E1_img_000000011_Default_000.tif';'180608_E1_img_000000014_Default_000.tif';'180608_E1_img_000000013_Default_000.tif';'180608_E1_img_000000016_Default_000.tif';'180608_E1_img_000000015_Default_000.tif';'180608_E1_img_000000018_Default_000.tif';'180608_E1_img_000000017_Default_000.tif'};
tags={'Green';'Red'};
finalfile=('180608_E1'); 
end
%% input control
    if size (before{1},1) ~= size (after{1},1)
     disp ('the number of images in the before and after conditions do not match!')
        return
    end

%% data import
    for lp=1:size(before,1)
        data.raw.before{lp} = imread(before{lp});
        data.raw.before{lp}=double(data.raw.before{lp});       
    end
    for a=1:size(after,1)
        data.raw.after{a} = imread(after{a});
        data.raw.after{a}=double(data.raw.after{a});      
    end
    


data.labels = tags; 
data.inputs.before = before; 
data.inputs.after = after;


%% Selection of region of interest for each neuron
first_green=horzcat(data.raw.before{1,1}(:));
first_red=horzcat(data.raw.before{1,2}(:));
first_green_max=max(first_green);
first_green_min=min(first_green);
first_red_max=max(first_red);
first_red_min=min(first_red);
brightest_red=horzcat(data.raw.after{1,(end-2)}(:));
brightest_red_max=max(brightest_red);


brighthening_coeff=255./first_green_max;
stack_proj_norm=imread(before{1,1});

if (nargin>1)
interval_spike_counts=input('What is the number of the spikes in each interval? (e.g. 250; !!if standard experiment enter 0) = ');
n=input('number of neurons to be examined =');
directory_to_save=input('Enter the name of the folder you wish to save the current file to (as a character variable with single quotations.) e.g. "Medium Sensitivity NB w patch"= ');
memory=input('Would you like to get the ROI values from the memory? (yes=1, no=0)=');

if n==1 & memory==1
    
    load('ROI_array_180608_E1_default.mat');
elseif n~=1 & memory==1
    oldFolder = cd('JSD'); 
    cd
    s1='JSD values';
    s2='.mat';
    filename=strcat(s1,finalfile,s2);
    matfileobject=matfile(filename);
    ROI_array=matfileobject.ROI_array;
    cd(oldFolder);
elseif n~=1 & memory==0
% first neuron is always the recorded neuron and the rest should be chosen
% from closest non recorded to furthest non recorded neuron+ 2 noise
% regions
for i=1:n
ROI = user_hand_drawn_ROI (stack_proj_norm.*brighthening_coeff);
ROI_array{i}=ROI;
end
end

end


%% plot data in time 

%recorded neuron roi extraction

    
    data.raw.beforerecorded{1,1}=data.raw.before{1,1}(ROI_array{1});
    data.raw.beforerecorded{1,2}=data.raw.before{1,2}(ROI_array{1});
    for i=1:size(after,1)
    data.raw.afterrecorded{1,i}=data.raw.after{1,i}(ROI_array{1});
    end
%first non recorded neuron (first post synaptic neuron) roi extraction
    data.raw.beforenonrecorded{1,1}=data.raw.before{1,1}(ROI_array{2});
    data.raw.beforenonrecorded{1,2}=data.raw.before{1,2}(ROI_array{2});

    for i=1:size(after,1)
    data.raw.afternonrecorded{1,i}=data.raw.after{1,i}(ROI_array{2});
    end
%second non recorded neuron (second post synaptic neuron) roi extraction
    data.raw.beforenonrecorded2{1,1}=data.raw.before{1,1}(ROI_array{3});
    data.raw.beforenonrecorded2{1,2}=data.raw.before{1,2}(ROI_array{3});
    
    for i=1:size(after,1)
    data.raw.afternonrecorded2{1,i}=data.raw.after{1,i}(ROI_array{3});
    end
%third non recorded neuron (third post synaptic neuron) roi extraction
    data.raw.beforenonrecorded3{1,1}=data.raw.before{1,1}(ROI_array{4});
    data.raw.beforenonrecorded3{1,2}=data.raw.before{1,2}(ROI_array{4});  
    for i=1:size(after,1)
    data.raw.afternonrecorded3{1,i}=data.raw.after{1,i}(ROI_array{4});
    end
     
%         
% %fourth non recorded neuron (fourth post synaptic neuron)roi extraction
%     
%     data.raw.beforenonrecorded4{1,1}=data.raw.before{1,1}(ROI_array{5});
%     data.raw.beforenonrecorded4{1,2}=data.raw.before{1,2}(ROI_array{5});
%     
%     for i=1:size(after,1)
%     data.raw.afternonrecorded4{1,i}=data.raw.after{1,i}(ROI_array{5});   
%     end
% 
% 
% %noise_1 ROI extraction
%     data.raw.beforenoise_1{1,1}=data.raw.before{1,1}(ROI_array{6});  
%     data.raw.beforenoise_1{1,2}=data.raw.before{1,2}(ROI_array{6});
%     
%     for i=1:size(after,1)
%     data.raw.afternoise_1{1,i}=data.raw.after{1,i}(ROI_array{6});
%     end
% %noise_2 ROI extraction
%     data.raw.beforenoise_2{1,1}=data.raw.before{1,1}(ROI_array{7});
%     data.raw.beforenoise_2{1,2}=data.raw.before{1,2}(ROI_array{7});
%     for i=1:size(after,1)
%     data.raw.afternoise_2{1,i}=data.raw.after{1,i}(ROI_array{7});
%     end


%% removing noise by exatracting the min value and normalizing with the max value of the highest pixel at the recorded neurons ROI.
%min value is taken as the general min
before_all=horzcat(data.raw.before{:});
before_all_horz=horzcat(before_all(:));
after_all=horzcat(data.raw.after{:});
after_all_horz=horzcat(after_all(:));
all_vert=[after_all_horz;before_all_horz];
ultimate_max=max(all_vert);
ultimate_min=min(all_vert);


% recorded neuron ROI red channel ultimate mx and min values to check
% before normalization
k=1;
for i=2:2:size(after,1)
    red_horizontal{k}=horzcat(data.raw.afterrecorded{1,i});
    R(1,k)=mean(red_horizontal{k});
    k=k+1;
end
red_all=vertcat(red_horizontal{:});
red_all_ROI_max=max(red_all);
red_all_ROI_min=min(red_all);
Rmax=max(R);
idx_Rmax=(R==Rmax);
Rmax_index=find(idx_Rmax==1);


%max pixel is taken from the recorded neuron ROI first image gren channel (in theory that is the brightest image)

y=horzcat(data.raw.beforerecorded{1,1});
first_green_ROI_max_pixel=max(y); %compare max pixel of first green with red ROI max pixel if red is bigger than green take red
max_normalize=first_green_ROI_max_pixel;
if (red_all_ROI_max>first_green_ROI_max_pixel)
    max_normalize=red_all_ROI_max;
end
    
% get the ratio of the ultimate max and first image recorded neuron ROI green
% channel max for further analysis on initial value effect
green_ROI_max_ultimate_max_ratio=first_green_ROI_max_pixel/ultimate_max;
green_ROI_initial_value_mean=mean(y);

%pixel to use in normalization

normalization_pixel=max_normalize-ultimate_min;


%% Normalization

% recorded neuron
data.normal.beforerecorded{1,1}=(data.raw.beforerecorded{1,1}-(ultimate_min))./normalization_pixel;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afterrecorded{1,d}=(data.raw.afterrecorded{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

data.normal.beforerecorded{1,2}=(data.raw.beforerecorded{1,2}-(ultimate_min))./normalization_pixel;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afterrecorded{1,d}=(data.raw.afterrecorded{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

%Normalizing the noise is different, depending on which channel is examined
%the max pixel will be taken from that channel normazlized by that.

% 
% % noise_1
% data.normal.beforenoise_1{1,1}=(data.raw.beforenoise_1{1,1}-(ultimate_min))./normalization_pixel;
% 
% d=1;
% m=size(after,1)/2;
% for l=1:m
%     data.normal.afternoise_1{1,d}=(data.raw.afternoise_1{1,d}-(ultimate_min))./normalization_pixel;
%     d=d+2;
% end
% 
% data.normal.beforenoise_1{1,2}=(data.raw.beforenoise_1{1,2}-(ultimate_min))./normalization_pixel;
% 
% d=2;
% m=size(after,1)/2;
% for l=1:m
%     data.normal.afternoise_1{1,d}=(data.raw.afternoise_1{1,d}-(ultimate_min))./normalization_pixel;
%     d=d+2;
% end
% 
% % noise_2
% data.normal.beforenoise_2{1,1}=(data.raw.beforenoise_2{1,1}-(ultimate_min))./normalization_pixel;
% 
% d=1;
% m=size(after,1)/2;
% for l=1:m
%     data.normal.afternoise_2{1,d}=(data.raw.afternoise_2{1,d}-(ultimate_min))./normalization_pixel;
%     d=d+2;
% end
% 
% data.normal.beforenoise_2{1,2}=(data.raw.beforenoise_2{1,2}-(ultimate_min))./normalization_pixel;
% 
% d=2;
% m=size(after,1)/2;
% for l=1:m
%     data.normal.afternoise_2{1,d}=(data.raw.afternoise_2{1,d}-(ultimate_min))./normalization_pixel;
%     d=d+2;
% end


%first non recorded neuron

data.normal.beforenonrecorded{1,1}=(data.raw.beforenonrecorded{1,1}-(ultimate_min))./normalization_pixel;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded{1,d}=(data.raw.afternonrecorded{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

data.normal.beforenonrecorded{1,2}=(data.raw.beforenonrecorded{1,2}-(ultimate_min))./normalization_pixel;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded{1,d}=(data.raw.afternonrecorded{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

%second non recorded neuron

data.normal.beforenonrecorded2{1,1}=(data.raw.beforenonrecorded2{1,1}-(ultimate_min))./normalization_pixel;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded2{1,d}=(data.raw.afternonrecorded2{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

data.normal.beforenonrecorded2{1,2}=(data.raw.beforenonrecorded2{1,2}-(ultimate_min))./normalization_pixel;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded2{1,d}=(data.raw.afternonrecorded2{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

%third non recorded neuron

data.normal.beforenonrecorded3{1,1}=(data.raw.beforenonrecorded3{1,1}-(ultimate_min))./normalization_pixel;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded3{1,d}=(data.raw.afternonrecorded3{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

data.normal.beforenonrecorded3{1,2}=(data.raw.beforenonrecorded3{1,2}-(ultimate_min))./normalization_pixel;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded3{1,d}=(data.raw.afternonrecorded3{1,d}-(ultimate_min))./normalization_pixel;
    d=d+2;
end

%fourth non recorded neuron
% 
% data.normal.beforenonrecorded4{1,1}=(data.raw.beforenonrecorded4{1,1}-(ultimate_min))./normalization_pixel;
% 
% d=1;
% m=size(after,1)/2;
% for l=1:m
%     data.normal.afternonrecorded4{1,d}=(data.raw.afternonrecorded4{1,d}-(ultimate_min))./normalization_pixel;
%     d=d+2;
% end
% 
% data.normal.beforenonrecorded4{1,2}=(data.raw.beforenonrecorded4{1,2}-(ultimate_min))./normalization_pixel;
% 
% d=2;
% m=size(after,1)/2;
% for l=1:m
%     data.normal.afternonrecorded4{1,d}=(data.raw.afternonrecorded4{1,d}-(ultimate_min))./normalization_pixel;
%     d=d+2;
% end


%% Calculating the difference between two conditions, through the change (increase or decrease) in percentage of each pixel between conditions.
%% pre - each condition comparison
% recorded neuron
% green channel recorded neuron
k=1;
for l=1:2:size(after,1)
    Green_pre_diff{k}=((abs(data.normal.beforerecorded{1,1}-data.normal.afterrecorded{1,l}))./data.normal.beforerecorded{1,1});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation. Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_pre_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_pre_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_pre_diff{i});
     Green_pre_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calculated
     % above.
     idx_inf = (Green_pre_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_pre_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_pre_diff{i}(~idx_inf));
        Green_pre_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_pre_diff_total(l)=(sum(Green_pre_diff{l}))./size(data.normal.beforerecorded{1,1},1);
end
% red channel recorded neuron
k=1;
for l=2:2:size(after,1)
    Red_pre_diff{k}=((abs(data.normal.beforerecorded{1,2}-data.normal.afterrecorded{1,l}))./data.normal.beforerecorded{1,2});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_pre_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_pre_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_pre_diff{i});
     Red_pre_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_pre_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_pre_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_pre_diff{i}(~idx_inf));
        Red_pre_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

m=size(after,1)/2;
for l=1:m
    Red_pre_diff_total(l)=(sum(Red_pre_diff{l}))./size(data.normal.beforerecorded{1,1},1);   
end

% conversion rate (R/G) recorded neuron
k=1;
for l=1:2:size(after,1)
    Conversion_pre_diff{k}=(abs((data.normal.beforerecorded{1,2}./data.normal.beforerecorded{1,1})-(data.normal.afterrecorded{1,l+1}./data.normal.afterrecorded{1,l})))./(data.normal.beforerecorded{1,2}./data.normal.beforerecorded{1,1});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_pre_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_pre_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_pre_diff{i});
     Conversion_pre_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_pre_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_pre_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_pre_diff{i}(~idx_inf));
        Conversion_pre_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_pre_diff_total(l)=(sum(Conversion_pre_diff{l}))./size(data.normal.beforerecorded{1,1},1);   
end
% first non recorded neuron
% green channel first nonrecorded neuron
k=1;
for l=1:2:size(after,1)
    Green_pre_N_diff{k}=((abs(data.normal.beforenonrecorded{1,1}-data.normal.afternonrecorded{1,l}))./data.normal.beforenonrecorded{1,1});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_pre_N_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_pre_N_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_pre_N_diff{i});
     Green_pre_N_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_pre_N_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_pre_N_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_pre_N_diff{i}(~idx_inf));
        Green_pre_N_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

m=size(after,1)/2;
for l=1:m
    Green_pre_N_diff_total(l)=(sum(Green_pre_N_diff{l}))./size(data.normal.beforenonrecorded{1,1},1);
end
% red channel first nonrecorded neuron
k=1;
for l=2:2:size(after,1)
    Red_pre_N_diff{k}=((abs(data.normal.beforenonrecorded{1,2}-data.normal.afternonrecorded{1,l}))./data.normal.beforenonrecorded{1,2});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_pre_N_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_pre_N_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_pre_N_diff{i});
     Red_pre_N_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_pre_N_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_pre_N_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_pre_N_diff{i}(~idx_inf));
        Red_pre_N_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_pre_N_diff_total(l)=(sum(Red_pre_N_diff{l}))./size(data.normal.beforenonrecorded{1,1},1);   
end

% conversion rate (R/G) first nonrecorded neuron
k=1;
for l=1:2:size(after,1)
    Conversion_pre_N_diff{k}=(abs((data.normal.beforenonrecorded{1,2}./data.normal.beforenonrecorded{1,1})-(data.normal.afternonrecorded{1,l+1}./data.normal.afternonrecorded{1,l})))./(data.normal.beforenonrecorded{1,2}./data.normal.beforenonrecorded{1,1});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_pre_N_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_pre_N_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_pre_N_diff{i});
     Conversion_pre_N_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_pre_N_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_pre_N_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_pre_N_diff{i}(~idx_inf));
        Conversion_pre_N_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_pre_N_diff_total(l)=(sum(Conversion_pre_N_diff{l}))./size(data.normal.beforenonrecorded{1,1},1);   
end

% second non recorded neuron
% green channel second nonrecorded neuron
k=1;
for l=1:2:size(after,1)
    Green_pre_N2_diff{k}=((abs(data.normal.beforenonrecorded2{1,1}-data.normal.afternonrecorded2{1,l}))./data.normal.beforenonrecorded2{1,1});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_pre_N2_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_pre_N2_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_pre_N2_diff{i});
     Green_pre_N2_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_pre_N2_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_pre_N2_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_pre_N2_diff{i}(~idx_inf));
        Green_pre_N2_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_pre_N2_diff_total(l)=(sum(Green_pre_N2_diff{l}))./size(data.normal.beforenonrecorded2{1,1},1);
end
% red channel second nonrecorded neuron
k=1;
for l=2:2:size(after,1)
    Red_pre_N2_diff{k}=((abs(data.normal.beforenonrecorded2{1,2}-data.normal.afternonrecorded2{1,l}))./data.normal.beforenonrecorded2{1,2});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_pre_N2_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_pre_N2_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_pre_N2_diff{i});
     Red_pre_N2_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_pre_N2_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_pre_N2_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_pre_N2_diff{i}(~idx_inf));
        Red_pre_N2_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_pre_N2_diff_total(l)=(sum(Red_pre_N2_diff{l}))./size(data.normal.beforenonrecorded2{1,1},1);   
end

% conversion rate (R/G) second nonrecorded neuron
k=1;
for l=1:2:size(after,1)
    Conversion_pre_N2_diff{k}=(abs((data.normal.beforenonrecorded2{1,2}./data.normal.beforenonrecorded2{1,1})-(data.normal.afternonrecorded2{1,l+1}./data.normal.afternonrecorded2{1,l})))./(data.normal.beforenonrecorded2{1,2}./data.normal.beforenonrecorded2{1,1});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_pre_N2_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_pre_N2_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_pre_N2_diff{i});
     Conversion_pre_N2_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_pre_N2_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_pre_N2_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_pre_N2_diff{i}(~idx_inf));
        Conversion_pre_N2_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;



m=size(after,1)/2;
for l=1:m
    Conversion_pre_N2_diff_total(l)=(sum(Conversion_pre_N2_diff{l}))./size(data.normal.beforenonrecorded2{1,1},1);   
end

% third non recorded neuron
% green channel third nonrecorded neuron
k=1;
for l=1:2:size(after,1)
    Green_pre_N3_diff{k}=((abs(data.normal.beforenonrecorded3{1,1}-data.normal.afternonrecorded3{1,l}))./data.normal.beforenonrecorded3{1,1});
    k=k+1;
end



% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_pre_N3_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_pre_N3_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_pre_N3_diff{i});
     Green_pre_N3_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_pre_N3_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_pre_N3_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_pre_N3_diff{i}(~idx_inf));
        Green_pre_N3_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_pre_N3_diff_total(l)=(sum(Green_pre_N3_diff{l}))./size(data.normal.beforenonrecorded3{1,1},1);
end
% red channel third nonrecorded neuron
k=1;
for l=2:2:size(after,1)
    Red_pre_N3_diff{k}=((abs(data.normal.beforenonrecorded3{1,2}-data.normal.afternonrecorded3{1,l}))./data.normal.beforenonrecorded3{1,2});
    k=k+1;
end



% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_pre_N3_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_pre_N3_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_pre_N3_diff{i});
     Red_pre_N3_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_pre_N3_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_pre_N3_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_pre_N3_diff{i}(~idx_inf));
        Red_pre_N3_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_pre_N3_diff_total(l)=(sum(Red_pre_N3_diff{l}))./size(data.normal.beforenonrecorded3{1,1},1);   
end

% conversion rate (R/G) third nonrecorded neuron
k=1;
for l=1:2:size(after,1)
    Conversion_pre_N3_diff{k}=(abs((data.normal.beforenonrecorded3{1,2}./data.normal.beforenonrecorded3{1,1})-(data.normal.afternonrecorded3{1,l+1}./data.normal.afternonrecorded3{1,l})))./(data.normal.beforenonrecorded3{1,2}./data.normal.beforenonrecorded3{1,1});
    k=k+1;
end



% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_pre_N3_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_pre_N3_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_pre_N3_diff{i});
     Conversion_pre_N3_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_pre_N3_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_pre_N3_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_pre_N3_diff{i}(~idx_inf));
        Conversion_pre_N3_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_pre_N3_diff_total(l)=(sum(Conversion_pre_N3_diff{l}))./size(data.normal.beforenonrecorded3{1,1},1);   
end
% fourth non recorded neuron
% green channel fourth nonrecorded neuron
% k=1;
% for l=1:2:size(after,1)
%     Green_pre_N4_diff{k}=((abs(data.normal.beforenonrecorded4{1,1}-data.normal.afternonrecorded4{1,l}))./data.normal.beforenonrecorded4{1,1});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Green_pre_N4_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Green_pre_N4_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Green_pre_N4_diff{i});
%      Green_pre_N4_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Green_pre_N4_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Green_pre_N4_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Green_pre_N4_diff{i}(~idx_inf));
%         Green_pre_N4_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Green_pre_N4_diff_total(l)=(sum(Green_pre_N4_diff{l}))./size(data.normal.beforenonrecorded4{1,1},1);
% end
% % red channel fourth nonrecorded neuron
% k=1;
% for l=2:2:size(after,1)
%     Red_pre_N4_diff{k}=((abs(data.normal.beforenonrecorded4{1,2}-data.normal.afternonrecorded4{1,l}))./data.normal.beforenonrecorded4{1,2});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Red_pre_N4_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Red_pre_N4_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Red_pre_N4_diff{i});
%      Red_pre_N4_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Red_pre_N4_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Red_pre_N4_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Red_pre_N4_diff{i}(~idx_inf));
%         Red_pre_N4_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Red_pre_N4_diff_total(l)=(sum(Red_pre_N4_diff{l}))./size(data.normal.beforenonrecorded4{1,1},1);   
% end
% 
% % conversion rate (R/G) fourth nonrecorded neuron
% k=1;
% for l=1:2:size(after,1)
%     Conversion_pre_N4_diff{k}=(abs((data.normal.beforenonrecorded4{1,2}./data.normal.beforenonrecorded4{1,1})-(data.normal.afternonrecorded4{1,l+1}./data.normal.afternonrecorded4{1,l})))./(data.normal.beforenonrecorded4{1,2}./data.normal.beforenonrecorded4{1,1});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Conversion_pre_N4_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Conversion_pre_N4_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Conversion_pre_N4_diff{i});
%      Conversion_pre_N4_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Conversion_pre_N4_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Conversion_pre_N4_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Conversion_pre_N4_diff{i}(~idx_inf));
%         Conversion_pre_N4_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Conversion_pre_N4_diff_total(l)=(sum(Conversion_pre_N4_diff{l}))./size(data.normal.beforenonrecorded4{1,1},1);   
% end
% 
% % noise 1
% % green channel noise 1
% k=1;
% for l=1:2:size(after,1)
%     Green_pre_noise1_diff{k}=((abs(data.normal.beforenoise_1{1,1}-data.normal.afternoise_1{1,l}))./data.normal.beforenoise_1{1,1});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Green_pre_noise1_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Green_pre_noise1_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Green_pre_noise1_diff{i});
%      Green_pre_noise1_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Green_pre_noise1_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Green_pre_noise1_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Green_pre_noise1_diff{i}(~idx_inf));
%         Green_pre_noise1_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Green_pre_noise1_diff_total(l)=(sum(Green_pre_noise1_diff{l}))./size(data.normal.beforenoise_1{1,1},1);
% end
% % red channel noise 1
% k=1;
% for l=2:2:size(after,1)
%     Red_pre_noise1_diff{k}=((abs(data.normal.beforenoise_1{1,2}-data.normal.afternoise_1{1,l}))./data.normal.beforenoise_1{1,2});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Red_pre_noise1_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Red_pre_noise1_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Red_pre_noise1_diff{i});
%      Red_pre_noise1_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Red_pre_noise1_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Red_pre_noise1_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Red_pre_noise1_diff{i}(~idx_inf));
%         Red_pre_noise1_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Red_pre_noise1_diff_total(l)=(sum(Red_pre_noise1_diff{l}))./size(data.normal.beforenoise_1{1,1},1);   
% end
% 
% % conversion rate (R/G) noise 1
% k=1;
% for l=1:2:size(after,1)
%     Conversion_pre_noise1_diff{k}=(abs((data.normal.beforenoise_1{1,2}./data.normal.beforenoise_1{1,1})-(data.normal.afternoise_1{1,l+1}./data.normal.afternoise_1{1,l})))./(data.normal.beforenoise_1{1,2}./data.normal.beforenoise_1{1,1});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Conversion_pre_noise1_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Conversion_pre_noise1_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Conversion_pre_noise1_diff{i});
%      Conversion_pre_noise1_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Conversion_pre_noise1_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Conversion_pre_noise1_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Conversion_pre_noise1_diff{i}(~idx_inf));
%         Conversion_pre_noise1_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Conversion_pre_noise1_diff_total(l)=(sum(Conversion_pre_noise1_diff{l}))./size(data.normal.beforenoise_1{1,1},1);   
% end
% 
% % noise 2
% % green channel noise 2
% k=1;
% for l=1:2:size(after,1)
%     Green_pre_noise2_diff{k}=((abs(data.normal.beforenoise_2{1,1}-data.normal.afternoise_2{1,l}))./data.normal.beforenoise_2{1,1});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Green_pre_noise2_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Green_pre_noise2_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Green_pre_noise2_diff{i});
%      Green_pre_noise2_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Green_pre_noise2_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Green_pre_noise2_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Green_pre_noise2_diff{i}(~idx_inf));
%         Green_pre_noise2_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Green_pre_noise2_diff_total(l)=(sum(Green_pre_noise2_diff{l}))./size(data.normal.beforenoise_2{1,1},1);
% end
% % red channel noise 2
% k=1;
% for l=2:2:size(after,1)
%     Red_pre_noise2_diff{k}=((abs(data.normal.beforenoise_2{1,2}-data.normal.afternoise_2{1,l}))./data.normal.beforenoise_2{1,2});
%     k=k+1;
% end
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Red_pre_noise2_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Red_pre_noise2_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Red_pre_noise2_diff{i});
%      Red_pre_noise2_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Red_pre_noise2_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Red_pre_noise2_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Red_pre_noise2_diff{i}(~idx_inf));
%         Red_pre_noise2_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Red_pre_noise2_diff_total(l)=(sum(Red_pre_noise2_diff{l}))./size(data.normal.beforenoise_2{1,1},1);   
% end
% 
% % conversion rate (R/G) noise 2
% k=1;
% for l=1:2:size(after,1)
%     Conversion_pre_noise2_diff{k}=(abs((data.normal.beforenoise_2{1,2}./data.normal.beforenoise_2{1,1})-(data.normal.afternoise_2{1,l+1}./data.normal.afternoise_2{1,l})))./(data.normal.beforenoise_2{1,2}./data.normal.beforenoise_2{1,1});
%     k=k+1;
% end
% 
% 
% 
% % Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% % which have Inf's as all elements
% all_values = horzcat(Conversion_pre_noise2_diff{:});
% idx_inf = (all_values == Inf);
% global_second_max = max(all_values(~idx_inf));
% 
% % second_max_alt = max(all_values(all_values ~=Inf));
% 
% % second_max_alt = max(all_values(all_values<max(all_values)));
% 
% for i = 1:size(Conversion_pre_noise2_diff,2)
% 
%      % Clean NaNs and assign 0's instead
%      idx_nan = isnan(Conversion_pre_noise2_diff{i});
%      Conversion_pre_noise2_diff{i}(idx_nan) = 0;
%         
%      % Clean Infs and assign the max value in the cell after the Inf, if
%      % all cell elements are Inf's assign the global max value calcukated
%      % above.
%      idx_inf = (Conversion_pre_noise2_diff{i} == Inf);
%      if nnz(idx_inf) == numel(idx_inf)
%          Conversion_pre_noise2_diff{i}(idx_inf) = global_second_max; 
%      else
%         local_second_max = max(Conversion_pre_noise2_diff{i}(~idx_inf));
%         Conversion_pre_noise2_diff{i}(idx_inf) = local_second_max; 
%      end
%          
% end
% clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% 
% 
% m=size(after,1)/2;
% for l=1:m
%     Conversion_pre_noise2_diff_total(l)=(sum(Conversion_pre_noise2_diff{l}))./size(data.normal.beforenoise_2{1,1},1);   
% end

%% consecutive condition comparison (pre-10spike)(10spike-20spike)...
% recorded neuron
% green channel recorded neuron

Green_consecutive_diff{1}=((abs(data.normal.beforerecorded{1,1}-data.normal.afterrecorded{1,1}))./data.normal.beforerecorded{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Green_consecutive_diff{k}=((abs(data.normal.afterrecorded{1,l+2}-data.normal.afterrecorded{1,l}))./data.normal.afterrecorded{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_consecutive_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_consecutive_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_consecutive_diff{i});
     Green_consecutive_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_consecutive_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_consecutive_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_consecutive_diff{i}(~idx_inf));
        Green_consecutive_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_consecutive_diff_total(l)=(sum(Green_consecutive_diff{l}))./size(data.normal.beforerecorded{1,1},1);
end
% red channel recorded neuron

Red_consecutive_diff{1}=((abs(data.normal.beforerecorded{1,2}-data.normal.afterrecorded{1,2}))./data.normal.beforerecorded{1,2});


k=2;
for l=2:2:(size(after,1)-2)
    Red_consecutive_diff{k}=((abs(data.normal.afterrecorded{1,l+2}-data.normal.afterrecorded{1,l}))./data.normal.afterrecorded{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_consecutive_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_consecutive_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_consecutive_diff{i});
     Red_consecutive_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_consecutive_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_consecutive_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_consecutive_diff{i}(~idx_inf));
        Red_consecutive_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_consecutive_diff_total(l)=(sum(Red_consecutive_diff{l}))./size(data.normal.beforerecorded{1,1},1);   
end

% conversion rate (R/G) recorded neuron
Conversion_consecutive_diff{1}=(abs((data.normal.beforerecorded{1,2}./data.normal.beforerecorded{1,1})-(data.normal.afterrecorded{1,2}./data.normal.afterrecorded{1,1})))./(data.normal.beforerecorded{1,2}./data.normal.beforerecorded{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Conversion_consecutive_diff{k}=(abs((data.normal.afterrecorded{1,l+1}./data.normal.afterrecorded{1,l})-(data.normal.afterrecorded{1,l+3}./data.normal.afterrecorded{1,l+2})))./(data.normal.afterrecorded{1,l+1}./data.normal.afterrecorded{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_consecutive_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_consecutive_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_consecutive_diff{i});
     Conversion_consecutive_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_consecutive_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_consecutive_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_consecutive_diff{i}(~idx_inf));
        Conversion_consecutive_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_consecutive_diff_total(l)=(sum(Conversion_consecutive_diff{l}))./size(data.normal.beforerecorded{1,1},1);   
end

% first nonrecorded neuron
% green channel first nonrecorded neuron

Green_consecutive_N_diff{1}=((abs(data.normal.beforenonrecorded{1,1}-data.normal.afternonrecorded{1,1}))./data.normal.beforenonrecorded{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Green_consecutive_N_diff{k}=((abs(data.normal.afternonrecorded{1,l+2}-data.normal.afternonrecorded{1,l}))./data.normal.afternonrecorded{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_consecutive_N_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_consecutive_N_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_consecutive_N_diff{i});
     Green_consecutive_N_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_consecutive_N_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_consecutive_N_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_consecutive_N_diff{i}(~idx_inf));
        Green_consecutive_N_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_consecutive_N_diff_total(l)=(sum(Green_consecutive_N_diff{l}))./size(data.normal.beforenonrecorded{1,1},1);
end
% red channel first nonrecorded neuron

Red_consecutive_N_diff{1}=((abs(data.normal.beforenonrecorded{1,2}-data.normal.afternonrecorded{1,2}))./data.normal.beforenonrecorded{1,2});


k=2;
for l=2:2:(size(after,1)-2)
    Red_consecutive_N_diff{k}=((abs(data.normal.afternonrecorded{1,l+2}-data.normal.afternonrecorded{1,l}))./data.normal.afternonrecorded{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_consecutive_N_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_consecutive_N_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_consecutive_N_diff{i});
     Red_consecutive_N_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_consecutive_N_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_consecutive_N_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_consecutive_N_diff{i}(~idx_inf));
        Red_consecutive_N_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_consecutive_N_diff_total(l)=(sum(Red_consecutive_N_diff{l}))./size(data.normal.beforenonrecorded{1,1},1);   
end

% conversion rate (R/G) first nonrecorded neuron
Conversion_consecutive_N_diff{1}=(abs((data.normal.beforenonrecorded{1,2}./data.normal.beforenonrecorded{1,1})-(data.normal.afternonrecorded{1,2}./data.normal.afternonrecorded{1,1})))./(data.normal.beforenonrecorded{1,2}./data.normal.beforenonrecorded{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Conversion_consecutive_N_diff{k}=(abs((data.normal.afternonrecorded{1,l+1}./data.normal.afternonrecorded{1,l})-(data.normal.afternonrecorded{1,l+3}./data.normal.afternonrecorded{1,l+2})))./(data.normal.afternonrecorded{1,l+1}./data.normal.afternonrecorded{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_consecutive_N_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_consecutive_N_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_consecutive_N_diff{i});
     Conversion_consecutive_N_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_consecutive_N_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_consecutive_N_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_consecutive_N_diff{i}(~idx_inf));
        Conversion_consecutive_N_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_consecutive_N_diff_total(l)=(sum(Conversion_consecutive_N_diff{l}))./size(data.normal.beforenonrecorded{1,1},1);   
end

% second nonrecorded neuron
% green channel second nonrecorded neuron

Green_consecutive_N2_diff{1}=((abs(data.normal.beforenonrecorded2{1,1}-data.normal.afternonrecorded2{1,1}))./data.normal.beforenonrecorded2{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Green_consecutive_N2_diff{k}=((abs(data.normal.afternonrecorded2{1,l+2}-data.normal.afternonrecorded2{1,l}))./data.normal.afternonrecorded2{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_consecutive_N2_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_consecutive_N2_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_consecutive_N2_diff{i});
     Green_consecutive_N2_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_consecutive_N2_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_consecutive_N2_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_consecutive_N2_diff{i}(~idx_inf));
        Green_consecutive_N2_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_consecutive_N2_diff_total(l)=(sum(Green_consecutive_N2_diff{l}))./size(data.normal.beforenonrecorded2{1,1},1);
end
% red channel second nonrecorded neuron

Red_consecutive_N2_diff{1}=((abs(data.normal.beforenonrecorded2{1,2}-data.normal.afternonrecorded2{1,2}))./data.normal.beforenonrecorded2{1,2});


k=2;
for l=2:2:(size(after,1)-2)
    Red_consecutive_N2_diff{k}=((abs(data.normal.afternonrecorded2{1,l+2}-data.normal.afternonrecorded2{1,l}))./data.normal.afternonrecorded2{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_consecutive_N2_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_consecutive_N2_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_consecutive_N2_diff{i});
     Red_consecutive_N2_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_consecutive_N2_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_consecutive_N2_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_consecutive_N2_diff{i}(~idx_inf));
        Red_consecutive_N2_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_consecutive_N2_diff_total(l)=(sum(Red_consecutive_N2_diff{l}))./size(data.normal.beforenonrecorded2{1,1},1);   
end

% conversion rate (R/G) second nonrecorded neuron
Conversion_consecutive_N2_diff{1}=(abs((data.normal.beforenonrecorded2{1,2}./data.normal.beforenonrecorded2{1,1})-(data.normal.afternonrecorded2{1,2}./data.normal.afternonrecorded2{1,1})))./(data.normal.beforenonrecorded2{1,2}./data.normal.beforenonrecorded2{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Conversion_consecutive_N2_diff{k}=(abs((data.normal.afternonrecorded2{1,l+1}./data.normal.afternonrecorded2{1,l})-(data.normal.afternonrecorded2{1,l+3}./data.normal.afternonrecorded2{1,l+2})))./(data.normal.afternonrecorded2{1,l+1}./data.normal.afternonrecorded2{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_consecutive_N2_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_consecutive_N2_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_consecutive_N2_diff{i});
     Conversion_consecutive_N2_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_consecutive_N2_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_consecutive_N2_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_consecutive_N2_diff{i}(~idx_inf));
        Conversion_consecutive_N2_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_consecutive_N2_diff_total(l)=(sum(Conversion_consecutive_N2_diff{l}))./size(data.normal.beforenonrecorded2{1,1},1);   
end

% third nonrecorded neuron
% green channel third nonrecorded neuron

Green_consecutive_N3_diff{1}=((abs(data.normal.beforenonrecorded3{1,1}-data.normal.afternonrecorded3{1,1}))./data.normal.beforenonrecorded3{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Green_consecutive_N3_diff{k}=((abs(data.normal.afternonrecorded3{1,l+2}-data.normal.afternonrecorded3{1,l}))./data.normal.afternonrecorded3{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_consecutive_N3_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_consecutive_N3_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_consecutive_N3_diff{i});
     Green_consecutive_N3_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_consecutive_N3_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_consecutive_N3_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_consecutive_N3_diff{i}(~idx_inf));
        Green_consecutive_N3_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Green_consecutive_N3_diff_total(l)=(sum(Green_consecutive_N3_diff{l}))./size(data.normal.beforenonrecorded3{1,1},1);
end
% red channel third nonrecorded neuron

Red_consecutive_N3_diff{1}=((abs(data.normal.beforenonrecorded3{1,2}-data.normal.afternonrecorded3{1,2}))./data.normal.beforenonrecorded3{1,2});


k=2;
for l=2:2:(size(after,1)-2)
    Red_consecutive_N3_diff{k}=((abs(data.normal.afternonrecorded3{1,l+2}-data.normal.afternonrecorded3{1,l}))./data.normal.afternonrecorded3{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_consecutive_N3_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_consecutive_N3_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_consecutive_N3_diff{i});
     Red_consecutive_N3_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_consecutive_N3_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_consecutive_N3_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_consecutive_N3_diff{i}(~idx_inf));
        Red_consecutive_N3_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Red_consecutive_N3_diff_total(l)=(sum(Red_consecutive_N3_diff{l}))./size(data.normal.beforenonrecorded3{1,1},1);   
end

% conversion rate (R/G) third nonrecorded neuron
Conversion_consecutive_N3_diff{1}=(abs((data.normal.beforenonrecorded3{1,2}./data.normal.beforenonrecorded3{1,1})-(data.normal.afternonrecorded3{1,2}./data.normal.afternonrecorded3{1,1})))./(data.normal.beforenonrecorded3{1,2}./data.normal.beforenonrecorded3{1,1});

k=2;
for l=1:2:(size(after,1)-2)
    Conversion_consecutive_N3_diff{k}=(abs((data.normal.afternonrecorded3{1,l+1}./data.normal.afternonrecorded3{1,l})-(data.normal.afternonrecorded3{1,l+3}./data.normal.afternonrecorded3{1,l+2})))./(data.normal.afternonrecorded3{1,l+1}./data.normal.afternonrecorded3{1,l});
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_consecutive_N3_diff{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_consecutive_N3_diff,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_consecutive_N3_diff{i});
     Conversion_consecutive_N3_diff{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_consecutive_N3_diff{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_consecutive_N3_diff{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_consecutive_N3_diff{i}(~idx_inf));
        Conversion_consecutive_N3_diff{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


m=size(after,1)/2;
for l=1:m
    Conversion_consecutive_N3_diff_total(l)=(sum(Conversion_consecutive_N3_diff{l}))./size(data.normal.beforenonrecorded3{1,1},1);   
end


%% all values in cell a cell array



all_diff_values{1}=Green_pre_diff_total;
all_diff_values{2}=Green_consecutive_diff_total;
all_diff_values{3}=Red_pre_diff_total;
all_diff_values{4}=Red_consecutive_diff_total;
all_diff_values{5}=Conversion_pre_diff_total;
all_diff_values{6}=Conversion_consecutive_diff_total;

all_diff_values{7}=Green_pre_N_diff_total;
all_diff_values{8}=Green_consecutive_N_diff_total;
all_diff_values{9}=Red_pre_N_diff_total;
all_diff_values{10}=Red_consecutive_N_diff_total;
all_diff_values{11}=Conversion_pre_N_diff_total;
all_diff_values{12}=Conversion_consecutive_N_diff_total;


all_diff_values{13}=Green_pre_N2_diff_total;
all_diff_values{14}=Green_consecutive_N2_diff_total;
all_diff_values{15}=Red_pre_N2_diff_total; 
all_diff_values{16}=Red_consecutive_N2_diff_total;
all_diff_values{17}=Conversion_pre_N2_diff_total;
all_diff_values{18}=Conversion_consecutive_N2_diff_total;

all_diff_values{19}=Green_pre_N3_diff_total;
all_diff_values{20}=Green_consecutive_N3_diff_total;
all_diff_values{21}=Red_pre_N3_diff_total;
all_diff_values{22}=Red_consecutive_N3_diff_total;
all_diff_values{23}=Conversion_pre_N3_diff_total;
all_diff_values{24}=Conversion_consecutive_N3_diff_total;

% all_diff_values{25}=[Red_consecutive_diff_total(1),Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2), Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2)+Red_consecutive_diff_total(3),Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2)+Red_consecutive_diff_total(3)+Red_consecutive_diff_total(4), Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2)+Red_consecutive_diff_total(3)+Red_consecutive_diff_total(4)+Red_consecutive_diff_total(5),Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2)+Red_consecutive_diff_total(3)+Red_consecutive_diff_total(4)+Red_consecutive_diff_total(5)+Red_consecutive_diff_total(6),Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2)+Red_consecutive_diff_total(3)+Red_consecutive_diff_total(4)+Red_consecutive_diff_total(5)+Red_consecutive_diff_total(6)+Red_consecutive_diff_total(7),Red_consecutive_diff_total(1)+Red_consecutive_diff_total(2)+Red_consecutive_diff_total(3)+Red_consecutive_diff_total(4)+Red_consecutive_diff_total(5)+Red_consecutive_diff_total(6)+Red_consecutive_diff_total(7)+Red_consecutive_diff_total(8)];
% all_diff_values{26}=[Red_consecutive_N_diff_total(1),Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2), Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2)+Red_consecutive_N_diff_total(3),Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2)+Red_consecutive_N_diff_total(3)+Red_consecutive_N_diff_total(4), Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2)+Red_consecutive_N_diff_total(3)+Red_consecutive_N_diff_total(4)+Red_consecutive_N_diff_total(5),Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2)+Red_consecutive_N_diff_total(3)+Red_consecutive_N_diff_total(4)+Red_consecutive_N_diff_total(5)+Red_consecutive_N_diff_total(6),Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2)+Red_consecutive_N_diff_total(3)+Red_consecutive_N_diff_total(4)+Red_consecutive_N_diff_total(5)+Red_consecutive_N_diff_total(6)+Red_consecutive_N_diff_total(7),Red_consecutive_N_diff_total(1)+Red_consecutive_N_diff_total(2)+Red_consecutive_N_diff_total(3)+Red_consecutive_N_diff_total(4)+Red_consecutive_N_diff_total(5)+Red_consecutive_N_diff_total(6)+Red_consecutive_N_diff_total(7)+Red_consecutive_N_diff_total(8)];
% all_diff_values{27}=[Red_consecutive_N2_diff_total(1),Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2), Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2)+Red_consecutive_N2_diff_total(3),Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2)+Red_consecutive_N2_diff_total(3)+Red_consecutive_N2_diff_total(4), Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2)+Red_consecutive_N2_diff_total(3)+Red_consecutive_N2_diff_total(4)+Red_consecutive_N2_diff_total(5),Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2)+Red_consecutive_N2_diff_total(3)+Red_consecutive_N2_diff_total(4)+Red_consecutive_N2_diff_total(5)+Red_consecutive_N2_diff_total(6),Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2)+Red_consecutive_N2_diff_total(3)+Red_consecutive_N2_diff_total(4)+Red_consecutive_N2_diff_total(5)+Red_consecutive_N2_diff_total(6)+Red_consecutive_N2_diff_total(7),Red_consecutive_N2_diff_total(1)+Red_consecutive_N2_diff_total(2)+Red_consecutive_N2_diff_total(3)+Red_consecutive_N2_diff_total(4)+Red_consecutive_N2_diff_total(5)+Red_consecutive_N2_diff_total(6)+Red_consecutive_N2_diff_total(7)+Red_consecutive_N2_diff_total(8)];
% all_diff_values{28}=[Red_consecutive_N3_diff_total(1),Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2), Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2)+Red_consecutive_N3_diff_total(3),Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2)+Red_consecutive_N3_diff_total(3)+Red_consecutive_N3_diff_total(4), Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2)+Red_consecutive_N3_diff_total(3)+Red_consecutive_N3_diff_total(4)+Red_consecutive_N3_diff_total(5),Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2)+Red_consecutive_N3_diff_total(3)+Red_consecutive_N3_diff_total(4)+Red_consecutive_N3_diff_total(5)+Red_consecutive_N3_diff_total(6),Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2)+Red_consecutive_N3_diff_total(3)+Red_consecutive_N3_diff_total(4)+Red_consecutive_N3_diff_total(5)+Red_consecutive_N3_diff_total(6)+Red_consecutive_N3_diff_total(7),Red_consecutive_N3_diff_total(1)+Red_consecutive_N3_diff_total(2)+Red_consecutive_N3_diff_total(3)+Red_consecutive_N3_diff_total(4)+Red_consecutive_N3_diff_total(5)+Red_consecutive_N3_diff_total(6)+Red_consecutive_N3_diff_total(7)+Red_consecutive_N3_diff_total(8)];

%% curve fitting
 ft=('linearinterp');

%fitting model separately to each diff total found

s=[];
if interval_spike_counts==0
    s=[10 30 80 180 430 880 1530];
else 
    i=0;
    for k=1:1:size(after,1)/2
    i=interval_spike_counts+i;
    s(k)=i;
    end
end
filename=strcat(directory_to_save,finalfile);
filename2=strcat(filename,'raw data multiplied by 100');
%prepares data to fit the curve, arranges data horizontally etc.
f1=figure; % raw data figure
f2=figure; % raw data multiplied_100 figure
for i=1:size(all_diff_values,2)  
% for l=1:7
% [P{i}{l}]=polyfit([s(l) s(l+1)],all_diff_values{i}(1,(l:l+1)),1);
% polynomf{l}=polyval(P{i}{l},[s(l) s(l+1)]);
% 
% %plot([s(l) s(l+1)],polynomf{l},'b');
% end
[xData, yData] = prepareCurveData( s, all_diff_values{i});
[fitresult{i}, gof{i}] = fit( xData, yData, ft);

if mod(i,6)==1
    color='g'; 
elseif mod(i,6)==2
    color='c';
elseif mod(i,6)==3
    color='m';
elseif mod(i,6)==4
    color='r';
elseif mod(i,6)==5
    color='k';
elseif mod(i,6)==0
    color='b';
end

figure(f1);
title(filename);
plot( fitresult{i}, color, xData, yData);
xlabel('Total action potentials fired');
ylabel('Average unit difference in 1 unit compared to the initial condition of the neuron')
hold on


% raw data multiplied with 100 to have the results in percentage
figure(f2);
title(filename2);
multiplied_100_yData=100.*(yData);
[fitresult_100{i}, gof_10{i}] = fit( xData, multiplied_100_yData, ft);
plot( fitresult_100{i}, color, xData, multiplied_100_yData);
xlabel('Total action potentials fired');
ylabel('Average difference in percentage compared to the initial condition of the neuron')
hold on
end
hold off
%% saving variables
% string=strcat(directory_to_save);
oldFolder = cd(directory_to_save); 



save([filename, '.mat']);
print(f1,filename,'-dsvg','-r600');
print(f2,filename2,'-dsvg','-r600');
%saveas(gcf,filename2,'png','-r300');
cd(oldFolder) ;

end
