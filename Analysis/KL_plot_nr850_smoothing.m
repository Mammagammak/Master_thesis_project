function KL_plot_nr850_smooothing(before,after,tags,finalfile)
%% this function analyzes the data using smoothing spline fitting model
close all
% this section is useful while debugging to prevent entering values over
% and over again
if (nargin < 1)
 
before = {'180608_E1_img_000000001_Default_000.tif';'180608_E1_img_000000002_Default_000.tif'};
after={'180608_E1_img_000000004_Default_000.tif';'180608_E1_img_000000003_Default_000.tif';'180608_E1_img_000000006_Default_000.tif';'180608_E1_img_000000005_Default_000.tif';'180608_E1_img_000000008_Default_000.tif';'180608_E1_img_000000007_Default_000.tif';'180608_E1_img_000000010_Default_000.tif';'180608_E1_img_000000009_Default_000.tif';'180608_E1_img_000000012_Default_000.tif';'180608_E1_img_000000011_Default_000.tif';'180608_E1_img_000000014_Default_000.tif';'180608_E1_img_000000013_Default_000.tif';'180608_E1_img_000000016_Default_000.tif';'180608_E1_img_000000015_Default_000.tif';'180608_E1_img_000000018_Default_000.tif';'180608_E1_img_000000017_Default_000.tif'}
tags={'Green';'Red'};
finalfile=('180608_E1');    

Non_recorded_neuron_count=4;
roi_recordedneuron_columns=[309 338];
roi_recordedneuron_rows=[386 418];
recorded_distance=0;


roi_nonrecordedneuron_columns=[260 296];
roi_nonrecordedneuron_rows=[349 394];
nonrecorded_distance=17.2;

roi_nonrecordedneuron2_columns=[472 502];
roi_nonrecordedneuron2_rows=[247 277];
nonrecorded2_distance=184.94;

roi_nonrecordedneuron3_columns=[67 111];
roi_nonrecordedneuron3_rows=[332 365];
nonrecorded3_distance=211.1;

roi_nonrecordedneuron4_columns=[477 509];
roi_nonrecordedneuron4_rows=[208 240];
nonrecorded4_distance=214.76;
   
    
end


%% Choosing ROI, the coordinates will be entered by the user individually for each neuron

if ~(nargin < 1)

Non_recorded_neuron_count=input('How many non recorded neurons were there in the experiment?=');

roi_recordedneuron_columns=input('enter the columns of the region of the recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_recordedneuron_rows=input('enter the rows of the region of the recorded neuron to analyze in an array. Exp:[y0 y1] = ');



roi_nonrecordedneuron_columns=input('enter the columns of the region of non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron_rows=input('enter the rows of the region of non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded_distance=input('enter the distance of the first non recorded neuron from the recorded neuron=');


roi_nonrecordedneuron2_columns=input('enter the columns of the region of the second non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron2_rows=input('enter the rows of the region of the second non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded2_distance=input('enter the distance of the second non recorded neuron from the recorded neuron=');


roi_nonrecordedneuron3_columns=input('enter the columns of the region of the third non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron3_rows=input('enter the rows of the region of the third non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded3_distance=input('enter the distance of the third non recorded neuron from the recorded neuron=');


roi_nonrecordedneuron4_columns=input('enter the columns of the region of the fourth non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron4_rows=input('enter the rows of the region of the fourth non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded4_distance=input('enter the distance of the fourth non recorded neuron from the recorded neuron=');

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



%% plot data in time 

%recorded neuron roi extraction

    
    data.raw.beforerecorded{1,1}=data.raw.before{1,1}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    data.raw.beforerecorded{1,2}=data.raw.before{1,2}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
  
    for i=1:size(after,1)
    data.raw.afterrecorded{1,i}=data.raw.after{1,i}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    end
%first non recorded neuron (first post synaptic neuron) roi extraction
    data.raw.beforenonrecorded{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    data.raw.beforenonrecorded{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));

    for i=1:size(after,1)
    data.raw.afternonrecorded{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    end
%second non recorded neuron (second post synaptic neuron) roi extraction
    data.raw.beforenonrecorded2{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    data.raw.beforenonrecorded2{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
  
    for i=1:size(after,1)
    data.raw.afternonrecorded2{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    end
%third non recorded neuron (third post synaptic neuron) roi extraction
    data.raw.beforenonrecorded3{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    data.raw.beforenonrecorded3{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
  
    for i=1:size(after,1)
    data.raw.afternonrecorded3{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    end
     
        
%fourth non recorded neuron (fourth post synaptic neuron)roi extraction
    
    data.raw.beforenonrecorded4{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
    data.raw.beforenonrecorded4{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
 
    for i=1:size(after,1)
    data.raw.afternonrecorded4{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
    end

    



    
%% For each ROI, max pixel value among the whole time series is found. The whole values are divided by this max value for the channel normalization.

% recorded neuron 

G(1,1)=mean2(data.raw.beforerecorded{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n
    G(1,b)=mean2(data.raw.afterrecorded{1,d}); 
     d=d+2;
     b=b+1;
    end
    Gmax=max(G);

    
R(1,1)=mean2(data.raw.beforerecorded{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m
    R(1,e)=mean2(data.raw.afterrecorded{1,f});
    f=f+2;
    e=e+1;
    end
    Rmax=max(R);


% first non recorded neuron 

G_N(1,1)=mean2(data.raw.beforenonrecorded{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n
    G_N(1,b)=mean2(data.raw.afternonrecorded{1,d}); 
     d=d+2;
     b=b+1;
    end
    G_Nmax=max(G_N);

    
R_N(1,1)=mean2(data.raw.beforenonrecorded{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m
    R_N(1,e)=mean2(data.raw.afternonrecorded{1,f});
    f=f+2;
    e=e+1;
    end
    R_Nmax=max(R_N);



% second non recorded neuron 

G_N2(1,1)=mean2(data.raw.beforenonrecorded2{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n
    G_N2(1,b)=mean2(data.raw.afternonrecorded2{1,d}); 
     d=d+2;
     b=b+1;
    end
    G_N2max=max(G_N2);


    
R_N2(1,1)=mean2(data.raw.beforenonrecorded2{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m
    R_N2(1,e)=mean2(data.raw.afternonrecorded2{1,f});
    f=f+2;
    e=e+1;
    end
    R_N2max=max(R_N2);


% third non recorded neuron 

G_N3(1,1)=mean2(data.raw.beforenonrecorded3{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n
    G_N3(1,b)=mean2(data.raw.afternonrecorded3{1,d}); 
     d=d+2;
     b=b+1;
    end
    G_N3max=max(G_N3);


    
R_N3(1,1)=mean2(data.raw.beforenonrecorded3{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m
    R_N3(1,e)=mean2(data.raw.afternonrecorded3{1,f});
    f=f+2;
    e=e+1;
    end
    R_N3max=max(R_N3);

% fourth non recorded neuron 

G_N4(1,1)=mean2(data.raw.beforenonrecorded4{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n
    G_N4(1,b)=mean2(data.raw.afternonrecorded4{1,d}); 
     d=d+2;
     b=b+1;
    end
    G_N4max=max(G_N4);


    
R_N4(1,1)=mean2(data.raw.beforenonrecorded4{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m
    R_N4(1,e)=mean2(data.raw.afternonrecorded4{1,f});
    f=f+2;
    e=e+1;
    end
    R_N4max=max(R_N4);

    
    
%% Normalizing the ROI's with the max value within each ROI

% recorded neuron
data.normal.beforerecorded{1,1}=data.raw.beforerecorded{1,1}./Gmax;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afterrecorded{1,d}=data.raw.afterrecorded{1,d}./Gmax;
    d=d+2;
end

data.normal.beforerecorded{1,2}=data.raw.beforerecorded{1,2}./Rmax;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afterrecorded{1,d}=data.raw.afterrecorded{1,d}./Rmax;
    d=d+2;
end

%first non recorded neuron

data.normal.beforenonrecorded{1,1}=data.raw.beforenonrecorded{1,1}./G_Nmax;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded{1,d}=data.raw.afternonrecorded{1,d}./G_Nmax;
    d=d+2;
end

data.normal.beforenonrecorded{1,2}=data.raw.beforenonrecorded{1,2}./R_Nmax;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded{1,d}=data.raw.afternonrecorded{1,d}./R_Nmax;
    d=d+2;
end

%second non recorded neuron

data.normal.beforenonrecorded2{1,1}=data.raw.beforenonrecorded2{1,1}./G_N2max;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded2{1,d}=data.raw.afternonrecorded2{1,d}./G_N2max;
    d=d+2;
end

data.normal.beforenonrecorded2{1,2}=data.raw.beforenonrecorded2{1,2}./R_N2max;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded2{1,d}=data.raw.afternonrecorded2{1,d}./R_N2max;
    d=d+2;
end

%third non recorded neuron

data.normal.beforenonrecorded3{1,1}=data.raw.beforenonrecorded3{1,1}./G_N3max;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded3{1,d}=data.raw.afternonrecorded3{1,d}./G_N3max;
    d=d+2;
end

data.normal.beforenonrecorded3{1,2}=data.raw.beforenonrecorded3{1,2}./R_N3max;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded3{1,d}=data.raw.afternonrecorded3{1,d}./R_N3max;
    d=d+2;
end

%foutrh non recorded neuron

data.normal.beforenonrecorded4{1,1}=data.raw.beforenonrecorded4{1,1}./G_N4max;

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded4{1,d}=data.raw.afternonrecorded4{1,d}./G_N4max;
    d=d+2;
end

data.normal.beforenonrecorded4{1,2}=data.raw.beforenonrecorded4{1,2}./R_N4max;

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded4{1,d}=data.raw.afternonrecorded4{1,d}./R_N4max;
    d=d+2;
end

%% converting the matrices into vectors before plotting the histogram
%recorded neuron
%green 
data.normal.beforerecordedvector{1,1}=data.normal.beforerecorded{1,1}(:);

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afterrecordedvector{1,d}=data.normal.afterrecorded{1,d}(:);
    d=d+2;
end

%red
data.normal.beforerecordedvector{1,2}=data.normal.beforerecorded{1,2}(:);

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afterrecordedvector{1,d}=data.normal.afterrecorded{1,d}(:);
    d=d+2;
end

%first non recorded neuron 
%green

data.normal.beforenonrecordedvector{1,1}=data.normal.beforenonrecorded{1,1}(:);

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecordedvector{1,d}=data.normal.afternonrecorded{1,d}(:);
    d=d+2;
end

%red
data.normal.beforenonrecordedvector{1,2}=data.normal.beforenonrecorded{1,2}(:);

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecordedvector{1,d}=data.normal.afternonrecorded{1,d}(:);
    d=d+2;
end

%second non recorded neuron 
%green

data.normal.beforenonrecorded2vector{1,1}=data.normal.beforenonrecorded2{1,1}(:);

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded2vector{1,d}=data.normal.afternonrecorded2{1,d}(:);
    d=d+2;
end

%red
data.normal.beforenonrecorded2vector{1,2}=data.normal.beforenonrecorded2{1,2}(:);

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded2vector{1,d}=data.normal.afternonrecorded2{1,d}(:);
    d=d+2;
end

%third non recorded neuron 
%green

data.normal.beforenonrecorded3vector{1,1}=data.normal.beforenonrecorded3{1,1}(:);

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded3vector{1,d}=data.normal.afternonrecorded3{1,d}(:);
    d=d+2;
end

%red
data.normal.beforenonrecorded3vector{1,2}=data.normal.beforenonrecorded3{1,2}(:);

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded3vector{1,d}=data.normal.afternonrecorded3{1,d}(:);
    d=d+2;
end

%fourth non recorded neuron 
%green

data.normal.beforenonrecorded4vector{1,1}=data.normal.beforenonrecorded4{1,1}(:);

d=1;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded4vector{1,d}=data.normal.afternonrecorded4{1,d}(:);
    d=d+2;
end

%red
data.normal.beforenonrecorded4vector{1,2}=data.normal.beforenonrecorded4{1,2}(:);

d=2;
m=size(after,1)/2;
for l=1:m
    data.normal.afternonrecorded4vector{1,d}=data.normal.afternonrecorded4{1,d}(:);
    d=d+2;
end

%% Calculating the Kernel Probability distribution for each condition for recorded neuron

%pts=[0.6, 0.669, 0.738, 0.807, 0.876, 0.945, 1.014, 1.083, 1.152, 1.221, 1.29];
pts=[0.2 0.28 0.36 0.44 0.52 0.6 0.68 0.76 0.84 0.92 1 1.08 1.16 1.24 1.32 1.4 1.48 1.56 1.64 1.72 1.8];    
% green channel
[fG{1},x]=ksdensity(data.normal.beforerecordedvector{1,1},pts);

d=1;
e=2;
m=size(after,1)/2;
for l=1:m
    [fG{e},x]=ksdensity(data.normal.afterrecordedvector{1,d},pts);
     d=d+2;
     e=e+1;
end

%red channel
[fR{1},x]=ksdensity(data.normal.beforerecordedvector{1,2},pts);

d=2;
e=2;
m=size(after,1)/2;
for l=1:m
    [fR{e},x]=ksdensity(data.normal.afterrecordedvector{1,d},pts);
     d=d+2;
     e=e+1;
end

% R/G conversion rate
% calculating the ratio and assigning them into new vectors
data.normal.recordedconversionratevector{1,1}=(data.normal.beforerecordedvector{1,2}./data.normal.beforerecordedvector{1,1});


m=size(after,1)/2;
d=1;
g=2;
for i=1:m
    data.normal.recordedconversionratevector{1,g}=(data.normal.afterrecordedvector{1,d+1}./data.normal.afterrecordedvector{1,(d)});
    d=d+2;
    g=g+1;
end

% clculationg kernel distribution
d=1;
e=1;
m=(size(after,1)/2)+1;
for l=1:m;
    [fC{e},x]=ksdensity(data.normal.recordedconversionratevector{1,d},pts);
     d=d+1;
     e=e+1;
end

% Green fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for i=1:m
        GreenpreKL{i}=fG{1,1}.*(log(fG{1,1}./fG{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(GreenpreKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(GreenpreKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(GreenpreKL{i});
     GreenpreKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (GreenpreKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         GreenpreKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(GreenpreKL{i}(~idx_inf));
        GreenpreKL{i}(idx_inf) = local_second_max; 
     end
         
end
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

a=1;
m=size(after,1)/2;

for l=1:m
    GreenpretotalKL(l)=sum(GreenpreKL{a});
    a=a+1;
end


% Red fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m;
        RedpreKL{i}=fR{1,1}.*(log(fR{1,1}./fR{1,k}));
    k=k+1;
    i=i+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(RedpreKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(RedpreKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(RedpreKL{i});
     RedpreKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (RedpreKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         RedpreKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(RedpreKL{i}(~idx_inf));
        RedpreKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m;
    RedpretotalKL(l)=sum(RedpreKL{a});
    a=a+1;
end

% Green fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 


m=size(after,1)/2;
k=1;
i=1;
for n=1:m;
    GreenconsecutiveKL{i}=fG{1,k}.*(log(fG{1,k}./fG{1,(k+1)}))
    k=k+1
    i=i+1
end

% Clean up of the NaN'a and Inf;s produced by the above calculation. Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(GreenconsecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(GreenconsecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(GreenconsecutiveKL{i});
     GreenconsecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (GreenconsecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         GreenconsecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(GreenconsecutiveKL{i}(~idx_inf));
        GreenconsecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m;
    GreenconsecutivetotalKL(l)=sum(GreenconsecutiveKL{l});
end


% Red fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
k=1;
i=1;
for n=1:m;
    RedconsecutiveKL{i}=fR{1,k}.*(log(fR{1,k}./fR{1,(k+1)}))
    k=k+1
    i=i+1
end

% Clean up of the NaN'a and Inf;s produced by the above calculation. Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(RedconsecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(RedconsecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(RedconsecutiveKL{i});
     RedconsecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (RedconsecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         RedconsecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(RedconsecutiveKL{i}(~idx_inf));
        RedconsecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m;
    RedconsecutivetotalKL(l)=sum(RedconsecutiveKL{l});
end

% R/G conversion rate KL divergence of each pre-x_spikes_condition


% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m;
        ConversionpreKL{n}=fC{1,1}.*(log(fC{1,1}./fC{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(ConversionpreKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(ConversionpreKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(ConversionpreKL{i});
     ConversionpreKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (ConversionpreKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         ConversionpreKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(ConversionpreKL{i}(~idx_inf));
        ConversionpreKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m;
    ConversionpretotalKL(l)=sum(ConversionpreKL{a});
    a=a+1;
end


% R/G conversion rate Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
k=1;
i=1;
for n=1:m;
    ConversionconsecutiveKL{i}=fC{1,k}.*(log(fC{1,k}./fC{1,(k+1)}))
    k=k+1
    i=i+1
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(ConversionconsecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(ConversionconsecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(ConversionconsecutiveKL{i});
     ConversionconsecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (ConversionconsecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         ConversionconsecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(ConversionconsecutiveKL{i}(~idx_inf));
        ConversionconsecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m;
    ConversionconsecutivetotalKL(l)=sum(ConversionconsecutiveKL{l});
end

%% Calculating the Kernel Probability distribution for each condition for first non recorded neuron


% green channel
[fG_N{1},x]=ksdensity(data.normal.beforenonrecordedvector{1,1},pts);

d=1;
e=2;
m=size(after,1)/2;
for l=1:m
    [fG_N{e},x]=ksdensity(data.normal.afternonrecordedvector{1,d},pts);
     d=d+2;
     e=e+1;
end

%red channel
[fR_N{1},x]=ksdensity(data.normal.beforenonrecordedvector{1,2},pts);

d=2;
e=2;
m=size(after,1)/2;
for l=1:m
    [fR_N{e},x]=ksdensity(data.normal.afternonrecordedvector{1,d},pts);
     d=d+2;
     e=e+1;
end

% R/G conversion rate
% calculating the ratio and assgning them into new vectors
data.normal.nonrecordedconversionratevector{1,1}=(data.normal.beforenonrecordedvector{1,2}./data.normal.beforenonrecordedvector{1,1});


m=size(after,1)/2;
d=1;
g=2;
for i=1:m
    data.normal.nonrecordedconversionratevector{1,g}=(data.normal.afternonrecordedvector{1,d+1}./data.normal.afternonrecordedvector{1,(d)});
    d=d+2;
    g=g+1;
end

% calculationg kernel distribution
d=1;
e=1;
m=(size(after,1)/2)+1;
for l=1:m
    [fC_N{e},x]=ksdensity(data.normal.nonrecordedconversionratevector{1,d},pts);
     d=d+1;
     e=e+1;
end

% Green fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m
        Green_NpreKL{i}=fG_N{1,1}.*(log(fG_N{1,1}./fG_N{1,k}));
    k=k+1;
    i=i+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_NpreKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_NpreKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_NpreKL{i});
     Green_NpreKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_NpreKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_NpreKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_NpreKL{i}(~idx_inf));
        Green_NpreKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution


m=size(after,1)/2;

for l=1:m
    Green_NpretotalKL(l)=sum(Green_NpreKL{l});
end




% Red fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m
        Red_NpreKL{i}=fR_N{1,1}.*(log(fR_N{1,1}./fR_N{1,k}));
    k=k+1;
    i=i+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_NpreKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_NpreKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_NpreKL{i});
     Red_NpreKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_NpreKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_NpreKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_NpreKL{i}(~idx_inf));
        Red_NpreKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Red_NpretotalKL(l)=sum(Red_NpreKL{a});
    a=a+1;
end

% Green fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 


m=size(after,1)/2;
k=1;
i=1;
for n=1:m
    Green_NconsecutiveKL{i}=fG_N{1,k}.*(log(fG_N{1,k}./fG_N{1,(k+1)}));
    k=k+1
    i=i+1
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_NconsecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_NconsecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_NconsecutiveKL{i});
     Green_NconsecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_NconsecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_NconsecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_NconsecutiveKL{i}(~idx_inf));
        Green_NconsecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Green_NconsecutivetotalKL(l)=sum(Green_NconsecutiveKL{l});
end


% Red fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
k=1;
i=1;
for n=1:m
    Red_NconsecutiveKL{i}=fR_N{1,k}.*(log(fR_N{1,k}./fR_N{1,(k+1)}));
    k=k+1
    i=i+1
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_NconsecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_NconsecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_NconsecutiveKL{i});
     Red_NconsecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_NconsecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_NconsecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_NconsecutiveKL{i}(~idx_inf));
        Red_NconsecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Red_NconsecutivetotalKL(l)=sum(Red_NconsecutiveKL{l});
end

% R/G conversion rate KL divergence of each pre-x_spikes_condition


% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Conversion_NpreKL{n}=fC_N{1,1}.*(log(fC_N{1,1}./fC_N{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_NpreKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_NpreKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_NpreKL{i});
     Conversion_NpreKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_NpreKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_NpreKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_NpreKL{i}(~idx_inf));
        Conversion_NpreKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Conversion_NpretotalKL(l)=sum(Conversion_NpreKL{a});
    a=a+1;
end


% R/G conversion rate Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
for n=1:m
    Conversion_NconsecutiveKL{n}=fC_N{1,n}.*(log(fC_N{1,n}./fC_N{1,(n+1)}));
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_NconsecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_NconsecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_NconsecutiveKL{i});
     Conversion_NconsecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_NconsecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_NconsecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_NconsecutiveKL{i}(~idx_inf));
        Conversion_NconsecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Conversion_NconsecutivetotalKL(l)=sum(Conversion_NconsecutiveKL{l});
end

%% Calculating the Kernel Probability distribution for each condition for second non recorded neuron



% green channel
[fG_N2{1},x]=ksdensity(data.normal.beforenonrecorded2vector{1,1},pts);

d=1;
e=2;
m=size(after,1)/2;
for l=1:m
    [fG_N2{e},x]=ksdensity(data.normal.afternonrecorded2vector{1,d},pts);
     d=d+2;
     e=e+1;
end

%red channel
[fR_N2{1},x]=ksdensity(data.normal.beforenonrecorded2vector{1,2},pts);

d=2;
e=2;
m=size(after,1)/2;
for l=1:m
    [fR_N2{e},x]=ksdensity(data.normal.afternonrecorded2vector{1,d},pts);
     d=d+2;
     e=e+1;
end

% R/G conversion rate
% calculating the ratio and assgning them into new vectors
data.normal.nonrecorded2conversionratevector{1,1}=(data.normal.beforenonrecorded2vector{1,2}./data.normal.beforenonrecorded2vector{1,1})


m=size(after,1)/2;
d=1;
g=2;
for i=1:m
    data.normal.nonrecorded2conversionratevector{1,g}=(data.normal.afternonrecorded2vector{1,d+1}./data.normal.afternonrecorded2vector{1,(d)})
    d=d+2;
    g=g+1;
end

% clculationg kernel distribution

m=(size(after,1)/2)+1;
for l=1:m
    [fC_N2{l},x]=ksdensity(data.normal.nonrecorded2conversionratevector{1,l},pts);
end

% Green fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m
        Green_N2preKL{i}=fG_N2{1,1}.*(log(fG_N2{1,1}./fG_N2{1,k}));
    k=k+1;
    i=i+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_N2preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_N2preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_N2preKL{i});
     Green_N2preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_N2preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_N2preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_N2preKL{i}(~idx_inf));
        Green_N2preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution


m=size(after,1)/2;

for l=1:m
    Green_N2pretotalKL(l)=sum(Green_N2preKL{l});
end




% Red fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m
        Red_N2preKL{i}=fR_N2{1,1}.*(log(fR_N2{1,1}./fR_N2{1,k}));
    k=k+1;
    i=i+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_N2preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_N2preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_N2preKL{i});
     Red_N2preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_N2preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_N2preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_N2preKL{i}(~idx_inf));
        Red_N2preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=size(after,1)/2;
for l=1:m
    Red_N2pretotalKL(l)=sum(Red_N2preKL{l});
end

% Green fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 


m=size(after,1)/2;
for n=1:m
    Green_N2consecutiveKL{n}=fG_N2{1,n}.*(log(fG_N2{1,n}./fG_N2{1,(n+1)}));
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_N2consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_N2consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_N2consecutiveKL{i});
     Green_N2consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_N2consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_N2consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_N2consecutiveKL{i}(~idx_inf));
        Green_N2consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Green_N2consecutivetotalKL(l)=sum(Green_N2consecutiveKL{l});
end


% Red fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
for n=1:m
    Red_N2consecutiveKL{n}=fR_N2{1,n}.*(log(fR_N2{1,n}./fR_N2{1,(n+1)}));
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_N2consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_N2consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_N2consecutiveKL{i});
     Red_N2consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_N2consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_N2consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_N2consecutiveKL{i}(~idx_inf));
        Red_N2consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Red_N2consecutivetotalKL(l)=sum(Red_N2consecutiveKL{l});
end

% R/G conversion rate KL divergence of each pre-x_spikes_condition


% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Conversion_N2preKL{n}=fC_N2{1,1}.*(log(fC_N2{1,1}./fC_N2{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_N2preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_N2preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_N2preKL{i});
     Conversion_N2preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_N2preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_N2preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_N2preKL{i}(~idx_inf));
        Conversion_N2preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Conversion_N2pretotalKL(l)=sum(Conversion_N2preKL{a});
    a=a+1;
end


% R/G conversion rate Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
for n=1:m
    Conversion_N2consecutiveKL{n}=fC_N2{1,n}.*(log(fC_N2{1,n}./fC_N2{1,(n+1)}));
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_N2consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_N2consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_N2consecutiveKL{i});
     Conversion_N2consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_N2consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Conversion_N2consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_N2consecutiveKL{i}(~idx_inf));
        Conversion_N2consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Conversion_N2consecutivetotalKL(l)=sum(Conversion_N2consecutiveKL{l});
end

%% Calculating the Kernel Probability distribution for each condition for third non recorded neuron



% green channel
[fG_N3{1},x]=ksdensity(data.normal.beforenonrecorded3vector{1,1},pts);

d=1;
e=2;
m=size(after,1)/2;
for l=1:m
    [fG_N3{e},x]=ksdensity(data.normal.afternonrecorded3vector{1,d},pts);
     d=d+2;
     e=e+1;
end

%red channel
[fR_N3{1},x]=ksdensity(data.normal.beforenonrecorded3vector{1,2},pts);

d=2;
e=2;
m=size(after,1)/2;
for l=1:m
    [fR_N3{e},x]=ksdensity(data.normal.afternonrecorded3vector{1,d},pts);
     d=d+2;
     e=e+1;
end

% R/G conversion rate
% calculating the ratio and assgning them into new vectors
data.normal.nonrecorded3conversionratevector{1,1}=(data.normal.beforenonrecorded3vector{1,2}./data.normal.beforenonrecorded3vector{1,1});


m=size(after,1)/2;
d=1;
g=2;
for i=1:m
    data.normal.nonrecorded3conversionratevector{1,g}=(data.normal.afternonrecorded3vector{1,d+1}./data.normal.afternonrecorded3vector{1,(d)});
    d=d+2;
    g=g+1;
end

% clculationg kernel distribution
d=1;
e=1;
m=(size(after,1)/2)+1;
for l=1:m
    [fC_N3{e},x]=ksdensity(data.normal.nonrecorded3conversionratevector{1,d},pts);
     d=d+1;
     e=e+1;
end

% Green fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Green_N3preKL{n}=fG_N3{1,1}.*(log(fG_N3{1,1}./fG_N3{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_N3preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_N3preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_N3preKL{i});
     Green_N3preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_N3preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Green_N3preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_N3preKL{i}(~idx_inf));
        Green_N3preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

a=1;
m=size(after,1)/2;

for l=1:m
    Green_N3pretotalKL(l)=sum(Green_N3preKL{a});
    a=a+1;
end




% Red fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Red_N3preKL{n}=fR_N3{1,1}.*(log(fR_N3{1,1}./fR_N3{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_N3preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_N3preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_N3preKL{i});
     Red_N3preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_N3preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
         Red_N3preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_N3preKL{i}(~idx_inf));
        Red_N3preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Red_N3pretotalKL(l)=sum(Red_N3preKL{a});
    a=a+1;
end

% Green fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 


m=size(after,1)/2;
k=1;
for n=1:m
    Green_N3consecutiveKL{n}=fG_N3{1,k}.*(log(fG_N3{1,k}./fG_N3{1,(k+1)}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_N3consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_N3consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_N3consecutiveKL{i});
     Green_N3consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_N3consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_N3consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_N3consecutiveKL{i}(~idx_inf));
        Green_N3consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Green_N3consecutivetotalKL(l)=sum(Green_N3consecutiveKL{l});
end


% Red fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
for n=1:m
    Red_N3consecutiveKL{n}=fR_N3{1,n}.*(log(fR_N3{1,n}./fR_N3{1,(n+1)}));
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_N3consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_N3consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_N3consecutiveKL{i});
     Red_N3consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_N3consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_N3consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_N3consecutiveKL{i}(~idx_inf));
        Red_N3consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;
% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Red_N3consecutivetotalKL(l)=sum(Red_N3consecutiveKL{l});
end

% R/G conversion rate KL divergence of each pre-x_spikes_condition


% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Conversion_N3preKL{n}=fC_N3{1,1}.*(log(fC_N3{1,1}./fC_N3{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_N3preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_N3preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_N3preKL{i});
     Conversion_N3preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_N3preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_N3preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_N3preKL{i}(~idx_inf));
        Conversion_N3preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Conversion_N3pretotalKL(l)=sum(Conversion_N3preKL{a});
    a=a+1;
end


% R/G conversion rate Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
for n=1:m
    Conversion_N3consecutiveKL{n}=fC_N3{1,n}.*(log(fC_N3{1,n}./fC_N3{1,(n+1)}));
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_N3consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_N3consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_N3consecutiveKL{i});
     Conversion_N3consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_N3consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_N3consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_N3consecutiveKL{i}(~idx_inf));
        Conversion_N3consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Conversion_N3consecutivetotalKL(l)=sum(Conversion_N3consecutiveKL{l});
end

%% Calculating the Kernel Probability distribution for each condition for fourth non recorded neuron



% green channel
[fG_N4{1},x]=ksdensity(data.normal.beforenonrecorded4vector{1,1},pts);

d=1;
e=2;
m=size(after,1)/2;
for l=1:m
    [fG_N4{e},x]=ksdensity(data.normal.afternonrecorded4vector{1,d},pts);
     d=d+2;
     e=e+1;
end

%red channel
[fR_N4{1},x]=ksdensity(data.normal.beforenonrecorded4vector{1,2},pts);

d=2;
e=2;
m=size(after,1)/2;
for l=1:m
    [fR_N4{e},x]=ksdensity(data.normal.afternonrecorded4vector{1,d},pts);
     d=d+2;
     e=e+1;
end

% R/G conversion rate
% calculating the ratio and assigning them into new vectors
data.normal.nonrecorded4conversionratevector{1,1}=(data.normal.beforenonrecorded4vector{1,2}./data.normal.beforenonrecorded4vector{1,1});


m=size(after,1)/2;
d=1;
g=2;
for i=1:m
    data.normal.nonrecorded4conversionratevector{1,g}=(data.normal.afternonrecorded4vector{1,d+1}./data.normal.afternonrecorded4vector{1,(d)});
    d=d+2;
    g=g+1;
end

% clculationg kernel distribution
d=1;
e=1;
m=(size(after,1)/2)+1;
for l=1:m
    [fC_N4{e},x]=ksdensity(data.normal.nonrecorded4conversionratevector{1,d},pts);
     d=d+1;
     e=e+1;
end

% Green fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Green_N4preKL{n}=fG_N4{1,1}.*(log(fG_N4{1,1}./fG_N4{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_N4preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_N4preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_N4preKL{i});
     Green_N4preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_N4preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_N4preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_N4preKL{i}(~idx_inf));
        Green_N4preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution


m=size(after,1)/2;

for l=1:m
    Green_N4pretotalKL(l)=sum(Green_N4preKL{l});
end




% Red fluorescence KL divergence of each pre-x_spikes_condition
% example: DKL(pre,10spikes), DKL(pre,20spikes)...

% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
i=1;
for n=1:m
        Red_N4preKL{i}=fR_N4{1,1}.*(log(fR_N4{1,1}./fR_N4{1,k}));
    k=k+1;
    i=i+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_N4preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_N4preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_N4preKL{i});
     Red_N4preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_N4preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_N4preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_N4preKL{i}(~idx_inf));
        Red_N4preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Red_N4pretotalKL(l)=sum(Red_N4preKL{a});
    a=a+1;
end

% Green fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 


m=size(after,1)/2;
k=1;
for n=1:m
    Green_N4consecutiveKL{n}=fG_N4{1,k}.*(log(fG_N4{1,k}./fG_N4{1,(k+1)}))
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_N4consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_N4consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_N4consecutiveKL{i});
     Green_N4consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_N4consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_N4consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_N4consecutiveKL{i}(~idx_inf));
        Green_N4consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Green_N4consecutivetotalKL(l)=sum(Green_N4consecutiveKL{l});
end


% Red fluorescence, Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
k=1;
for n=1:m
    Red_N4consecutiveKL{n}=fR_N4{1,k}.*(log(fR_N4{1,k}./fR_N4{1,(k+1)}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_N4consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_N4consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_N4consecutiveKL{i});
     Red_N4consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_N4consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_N4consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_N4consecutiveKL{i}(~idx_inf));
        Red_N4consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m
    Red_N4consecutivetotalKL(l)=sum(Red_N4consecutiveKL{l});
end

% R/G conversion rate KL divergence of each pre-x_spikes_condition


% this part will calculate the vectors for each couple compared condition 
m=size(after,1)/2;
k=2;
for n=1:m
        Conversion_N4preKL{n}=fC_N4{1,1}.*(log(fC_N4{1,1}./fC_N4{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_N4preKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_N4preKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_N4preKL{i});
     Conversion_N4preKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_N4preKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_N4preKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_N4preKL{i}(~idx_inf));
        Conversion_N4preKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=size(after,1)/2;
for l=1:m
    Conversion_N4pretotalKL(l)=sum(Conversion_N4preKL{l});
end


% R/G conversion rate Kl divergence of each condition with its previous condition
% example: DKL(pre,10spikes), DKL(10spikes,20spikes)...

% this part will calculate one vector for each couple compared condition 
m=size(after,1)/2;
k=1;

for n=1:m
    Conversion_N4consecutiveKL{n}=fC_N4{1,k}.*(log(fC_N4{1,k}./fC_N4{1,(k+1)}))
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_N4consecutiveKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_N4consecutiveKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_N4consecutiveKL{i});
     Conversion_N4consecutiveKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_N4consecutiveKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_N4consecutiveKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_N4consecutiveKL{i}(~idx_inf));
        Conversion_N4consecutiveKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
m=size(after,1)/2;
for l=1:m;
    Conversion_N4consecutivetotalKL(l)=sum(Conversion_N4consecutiveKL{l});
end

%% KL Divergence of recorded and non recorded comparisons
% recorded - first non recorded neuron each condition compared with its
% equal. exp: KL divergence of nonrecorded pre condition distribution and
% recorded pre condition distribution

%green corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Green_rec_N_correspondingKL{n}=fG{1,k}.*(log(fG{1,k}./fG_N{1,k}));
    k=k+1;
end

% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_rec_N_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_rec_N_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_rec_N_correspondingKL{i});
     Green_rec_N_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_rec_N_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_rec_N_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_rec_N_correspondingKL{i}(~idx_inf));
        Green_rec_N_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Green_rec_N_correspondingtotalKL(l)=sum(Green_rec_N_correspondingKL{l});
end

%red corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Red_rec_N_correspondingKL{n}=fR{1,k}.*(log(fR{1,k}./fR_N{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_rec_N_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_rec_N_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_rec_N_correspondingKL{i});
     Red_rec_N_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_rec_N_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_rec_N_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_rec_N_correspondingKL{i}(~idx_inf));
        Red_rec_N_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Red_rec_N_correspondingtotalKL(l)=sum(Red_rec_N_correspondingKL{l});
end

%conversion rate corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Conversion_rec_N_correspondingKL{n}=fC{1,k}.*(log(fC{1,k}./fC_N{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_rec_N_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_rec_N_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_rec_N_correspondingKL{i});
     Conversion_rec_N_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_rec_N_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_rec_N_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_rec_N_correspondingKL{i}(~idx_inf));
        Conversion_rec_N_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Conversion_rec_N_correspondingtotalKL(l)=sum(Conversion_rec_N_correspondingKL{l});
end


% recorded - second non recorded neuron each condition compared with its
% equal. exp: KL divergence of nonrecorded pre condition distribution and
% recorded pre condition distribution

%green corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Green_rec_N2_correspondingKL{n}=fG{1,k}.*(log(fG{1,k}./fG_N2{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_rec_N2_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_rec_N2_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_rec_N2_correspondingKL{i});
     Green_rec_N2_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_rec_N2_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_rec_N2_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_rec_N2_correspondingKL{i}(~idx_inf));
        Green_rec_N2_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Green_rec_N2_correspondingtotalKL(l)=sum(Green_rec_N2_correspondingKL{l});
end

%red corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Red_rec_N2_correspondingKL{n}=fR{1,k}.*(log(fR{1,k}./fR_N2{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_rec_N2_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_rec_N2_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_rec_N2_correspondingKL{i});
     Red_rec_N2_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_rec_N2_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_rec_N2_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_rec_N2_correspondingKL{i}(~idx_inf));
        Red_rec_N2_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Red_rec_N2_correspondingtotalKL(l)=sum(Red_rec_N2_correspondingKL{l});
end

%conversion rate corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Conversion_rec_N2_correspondingKL{n}=fC{1,k}.*(log(fC{1,k}./fC_N2{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_rec_N2_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_rec_N2_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_rec_N2_correspondingKL{i});
     Conversion_rec_N2_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_rec_N2_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_rec_N2_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_rec_N2_correspondingKL{i}(~idx_inf));
        Conversion_rec_N2_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Conversion_rec_N2_correspondingtotalKL(l)=sum(Conversion_rec_N2_correspondingKL{l});
end


% recorded - third non recorded neuron each condition compared with its
% equal. exp: KL divergence of nonrecorded pre condition distribution and
% recorded pre condition distribution

%green corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Green_rec_N3_correspondingKL{n}=fG{1,k}.*(log(fG{1,k}./fG_N3{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_rec_N3_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_rec_N3_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_rec_N3_correspondingKL{i});
     Green_rec_N3_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_rec_N3_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_rec_N3_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_rec_N3_correspondingKL{i}(~idx_inf));
        Green_rec_N3_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Green_rec_N3_correspondingtotalKL(l)=sum(Green_rec_N3_correspondingKL{l});
end

%red corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Red_rec_N3_correspondingKL{n}=fR{1,k}.*(log(fR{1,k}./fR_N3{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_rec_N3_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_rec_N3_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_rec_N3_correspondingKL{i});
     Red_rec_N3_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_rec_N3_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_rec_N3_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_rec_N3_correspondingKL{i}(~idx_inf));
        Red_rec_N3_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Red_rec_N3_correspondingtotalKL(l)=sum(Red_rec_N3_correspondingKL{l});
end

%conversion rate corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Conversion_rec_N3_correspondingKL{n}=fC{1,k}.*(log(fC{1,k}./fC_N3{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_rec_N3_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_rec_N3_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_rec_N3_correspondingKL{i});
     Conversion_rec_N3_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_rec_N3_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_rec_N3_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_rec_N3_correspondingKL{i}(~idx_inf));
        Conversion_rec_N3_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Conversion_rec_N3_correspondingtotalKL(l)=sum(Conversion_rec_N3_correspondingKL{l});
end


% recorded - fourth non recorded neuron each condition compared with its
% equal. exp: KL divergence of nonrecorded pre condition distribution and
% recorded pre condition distribution

%green corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Green_rec_N4_correspondingKL{n}=fG{1,k}.*(log(fG{1,k}./fG_N4{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Green_rec_N4_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Green_rec_N4_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Green_rec_N4_correspondingKL{i});
     Green_rec_N4_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Green_rec_N4_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Green_rec_N4_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Green_rec_N4_correspondingKL{i}(~idx_inf));
        Green_rec_N4_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Green_rec_N4_correspondingtotalKL(l)=sum(Green_rec_N4_correspondingKL{l});
end

%red corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Red_rec_N4_correspondingKL{n}=fR{1,k}.*(log(fR{1,k}./fR_N4{1,k}));
    k=k+1;
end


% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Red_rec_N4_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Red_rec_N4_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Red_rec_N4_correspondingKL{i});
     Red_rec_N4_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Red_rec_N4_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Red_rec_N4_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Red_rec_N4_correspondingKL{i}(~idx_inf));
        Red_rec_N4_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;


% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Red_rec_N4_correspondingtotalKL(l)=sum(Red_rec_N4_correspondingKL{l});
end

%conversion rate corresponding conditions


% this part will calculate the vectors for each couple compared condition 
m=(size(after,1)/2)+1;
k=1;
for n=1:m
        Conversion_rec_N4_correspondingKL{n}=fC{1,k}.*(log(fC{1,k}./fC_N4{1,k}));
    k=k+1;
end



% Clean up of the NaN'a and Inf;s produced by the above calculation.Find global maximum in order to assign to cells
% which have Inf's as all elements
all_values = horzcat(Conversion_rec_N4_correspondingKL{:});
idx_inf = (all_values == Inf);
global_second_max = max(all_values(~idx_inf));

% second_max_alt = max(all_values(all_values ~=Inf));

% second_max_alt = max(all_values(all_values<max(all_values)));

for i = 1:size(Conversion_rec_N4_correspondingKL,2)

     % Clean NaNs and assign 0's instead
     idx_nan = isnan(Conversion_rec_N4_correspondingKL{i});
     Conversion_rec_N4_correspondingKL{i}(idx_nan) = 0;
        
     % Clean Infs and assign the max value in the cell after the Inf, if
     % all cell elements are Inf's assign the global max value calcukated
     % above.
     idx_inf = (Conversion_rec_N4_correspondingKL{i} == Inf);
     if nnz(idx_inf) == numel(idx_inf)
        Conversion_rec_N4_correspondingKL{i}(idx_inf) = global_second_max; 
     else
        local_second_max = max(Conversion_rec_N4_correspondingKL{i}(~idx_inf));
        Conversion_rec_N4_correspondingKL{i}(idx_inf) = local_second_max; 
     end        
end
       
clearvars idx_inf idx_nan local_second_max global_second_max all_values;

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution

m=(size(after,1)/2)+1;
for l=1:m
    Conversion_rec_N4_correspondingtotalKL(l)=sum(Conversion_rec_N4_correspondingKL{l});
end
%% KL divergence calculated with subtracted Kernel distributions
% first non recorded neuron

% each neuron's each condition's distribution subtracted from its pre
% condition dstribution. exp: distributions subtract(pre-10spike, pre-20spike, pre-50spike)
% green, red and conversion rate
m=size(after,1)/2;
k=2;
for i=1:m
    fG_pre_difference{i}=fG{1,1}-fG{1,k};
    fG_Npre_difference{i}=fG_N{1,1}-fG_N{1,k};
    fG_N2pre_difference{i}=fG_N2{1,1}-fG_N2{1,k};
    fG_N3pre_difference{i}=fG_N3{1,1}-fG_N3{1,k};
    fG_N4pre_difference{i}=fG_N4{1,1}-fG_N4{1,k};
    fR_pre_difference{i}=fR{1,1}-fR{1,k};
    fR_Npre_difference{i}=fR_N{1,1}-fR_N{1,k};
    fR_N2pre_difference{i}=fR_N2{1,1}-fR_N2{1,k};
    fR_N3pre_difference{i}=fR_N3{1,1}-fR_N3{1,k};
    fR_N4pre_difference{i}=fR_N4{1,1}-fR_N4{1,k};
    fC_pre_difference{i}=fC{1,1}-fC{1,k};
    fC_Npre_difference{i}=fC_N{1,1}-fC_N{1,k};
    fC_N2pre_difference{i}=fC_N2{1,1}-fC_N2{1,k};
    fC_N3pre_difference{i}=fC_N3{1,1}-fC_N3{1,k};
    fC_N4pre_difference{i}=fC_N4{1,1}-fC_N4{1,k};
    k=k+1;
end
% each neuron's each condition's distribution subtracted from its previous
% condition dstribution. exp: distributions subtract(pre-10spike, 10spike-20spike, 20spike-50spike)
% green, red and conversion rate

m=size(after,1)/2;
k=1;
for i=1:m
    fG_consecutive_difference{i}=fG{1,k}-fG{1,(k+1)};
    fG_Nconsecutive_difference{i}=fG_N{1,k}-fG_N{1,(k+1)};
    fG_N2consecutive_difference{i}=fG_N2{1,k}-fG_N2{1,(k+1)};
    fG_N3consecutive_difference{i}=fG_N3{1,k}-fG_N3{1,(k+1)};
    fG_N4consecutive_difference{i}=fG_N4{1,k}-fG_N4{1,(k+1)};
    fR_consecutive_difference{i}=fR{1,k}-fR{1,(k+1)};
    fR_Nconsecutive_difference{i}=fR_N{1,k}-fR_N{1,(k+1)};
    fR_N2consecutive_difference{i}=fR_N2{1,k}-fR_N2{1,(k+1)};
    fR_N3consecutive_difference{i}=fR_N3{1,k}-fR_N3{1,(k+1)};
    fR_N4consecutive_difference{i}=fR_N4{1,k}-fR_N4{1,(k+1)};
    fC_consecutive_difference{i}=fC{1,k}-fC{1,(k+1)};
    fC_Nconsecutive_difference{i}=fC_N{1,k}-fC_N{1,(k+1)};
    fC_N2consecutive_difference{i}=fC_N2{1,k}-fC_N2{1,(k+1)};
    fC_N3consecutive_difference{i}=fC_N3{1,k}-fC_N3{1,(k+1)};
    fC_N4consecutive_difference{i}=fC_N4{1,k}-fC_N4{1,(k+1)};
    k=k+1;
end



% this part will calculate the vectors for each couple compared condition
% for pre-difference 
m=size(after,1)/2;
k=1;
i=1;
for n=1:m
        Green_rec_N_pre_differenceKL{n}=fG_pre_difference{1,k}.*(log(fG_pre_difference{1,k}./fG_Npre_difference{1,k}));
        Green_rec_N2_pre_differenceKL{n}=fG_pre_difference{1,k}.*(log(fG_pre_difference{1,k}./fG_N2pre_difference{1,k}));
        Green_rec_N3_pre_differenceKL{n}=fG_pre_difference{1,k}.*(log(fG_pre_difference{1,k}./fG_N3pre_difference{1,k}));
        Green_rec_N4_pre_differenceKL{n}=fG_pre_difference{1,k}.*(log(fG_pre_difference{1,k}./fG_N4pre_difference{1,k}));
        Red_rec_N_pre_differenceKL{n}=fR_pre_difference{1,k}.*(log(fR_pre_difference{1,k}./fR_Npre_difference{1,k}));
        Red_rec_N2_pre_differenceKL{n}=fR_pre_difference{1,k}.*(log(fR_pre_difference{1,k}./fR_N2pre_difference{1,k}));
        Red_rec_N3_pre_differenceKL{n}=fR_pre_difference{1,k}.*(log(fR_pre_difference{1,k}./fR_N3pre_difference{1,k}));
        Red_rec_N4_pre_differenceKL{n}=fR_pre_difference{1,k}.*(log(fR_pre_difference{1,k}./fR_N4pre_difference{1,k}));
        Conversion_rec_N_pre_differenceKL{n}=fC_pre_difference{1,k}.*(log(fC_pre_difference{1,k}./fC_Npre_difference{1,k}));
        Conversion_rec_N2_pre_differenceKL{n}=fC_pre_difference{1,k}.*(log(fC_pre_difference{1,k}./fC_N2pre_difference{1,k}));
        Conversion_rec_N3_pre_differenceKL{n}=fC_pre_difference{1,k}.*(log(fC_pre_difference{1,k}./fC_N3pre_difference{1,k}));
        Conversion_rec_N4_pre_differenceKL{n}=fC_pre_difference{1,k}.*(log(fC_pre_difference{1,k}./fC_N4pre_difference{1,k}));

        k=k+1;
end

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Green_rec_N_pre_differencetotalKL(l)=sum(Green_rec_N_pre_differenceKL{a});
    Green_rec_N2_pre_differencetotalKL(l)=sum(Green_rec_N2_pre_differenceKL{a});
    Green_rec_N3_pre_differencetotalKL(l)=sum(Green_rec_N3_pre_differenceKL{a});
    Green_rec_N4_pre_differencetotalKL(l)=sum(Green_rec_N4_pre_differenceKL{a});
    Red_rec_N_pre_differencetotalKL(l)=sum(Red_rec_N_pre_differenceKL{a});
    Red_rec_N2_pre_differencetotalKL(l)=sum(Red_rec_N2_pre_differenceKL{a});
    Red_rec_N3_pre_differencetotalKL(l)=sum(Red_rec_N3_pre_differenceKL{a});
    Red_rec_N4_pre_differencetotalKL(l)=sum(Red_rec_N4_pre_differenceKL{a});
    Conversion_rec_N_pre_differencetotalKL(l)=sum(Conversion_rec_N_pre_differenceKL{a});
    Conversion_rec_N2_pre_differencetotalKL(l)=sum(Conversion_rec_N2_pre_differenceKL{a});
    Conversion_rec_N3_pre_differencetotalKL(l)=sum(Conversion_rec_N3_pre_differenceKL{a});
    Conversion_rec_N4_pre_differencetotalKL(l)=sum(Conversion_rec_N4_pre_differenceKL{a});
   
    a=a+1;
end

% for consecutive condition pre-10spikes, 10-20spikes, 20-50spike....
m=size(after,1)/2;
k=1;
i=1;
for n=1:m
        Green_rec_N_consecutive_differenceKL{n}=fG_consecutive_difference{1,k}.*(log(fG_consecutive_difference{1,k}./fG_Nconsecutive_difference{1,k}));
        Green_rec_N2_consecutive_differenceKL{n}=fG_consecutive_difference{1,k}.*(log(fG_consecutive_difference{1,k}./fG_N2consecutive_difference{1,k}));
        Green_rec_N3_consecutive_differenceKL{n}=fG_consecutive_difference{1,k}.*(log(fG_consecutive_difference{1,k}./fG_N3consecutive_difference{1,k}));
        Green_rec_N4_consecutive_differenceKL{n}=fG_consecutive_difference{1,k}.*(log(fG_consecutive_difference{1,k}./fG_N4consecutive_difference{1,k}));
        Red_rec_N_consecutive_differenceKL{n}=fR_consecutive_difference{1,k}.*(log(fR_consecutive_difference{1,k}./fR_Nconsecutive_difference{1,k}));
        Red_rec_N2_consecutive_differenceKL{n}=fR_consecutive_difference{1,k}.*(log(fR_consecutive_difference{1,k}./fR_N2consecutive_difference{1,k}));
        Red_rec_N3_consecutive_differenceKL{n}=fR_consecutive_difference{1,k}.*(log(fR_consecutive_difference{1,k}./fR_N3consecutive_difference{1,k}));
        Red_rec_N4_consecutive_differenceKL{n}=fR_consecutive_difference{1,k}.*(log(fR_consecutive_difference{1,k}./fR_N4consecutive_difference{1,k}));
        Conversion_rec_N_consecutive_differenceKL{n}=fC_consecutive_difference{1,k}.*(log(fC_consecutive_difference{1,k}./fC_Nconsecutive_difference{1,k}));
        Conversion_rec_N2_consecutive_differenceKL{n}=fC_consecutive_difference{1,k}.*(log(fC_consecutive_difference{1,k}./fC_N2consecutive_difference{1,k}));
        Conversion_rec_N3_consecutive_differenceKL{n}=fC_consecutive_difference{1,k}.*(log(fC_consecutive_difference{1,k}./fC_N3consecutive_difference{1,k}));
        Conversion_rec_N4_consecutive_differenceKL{n}=fC_consecutive_difference{1,k}.*(log(fC_consecutive_difference{1,k}./fC_N4consecutive_difference{1,k}));

        k=k+1;
end

% this part will sum each individual vector within itself to have one final
% KL divergence value for each compared distribution
a=1;
m=size(after,1)/2;
for l=1:m
    Green_rec_N_consecutive_differencetotalKL(l)=sum(Green_rec_N_consecutive_differenceKL{a});
    Green_rec_N2_consecutive_differencetotalKL(l)=sum(Green_rec_N2_consecutive_differenceKL{a});
    Green_rec_N3_consecutive_differencetotalKL(l)=sum(Green_rec_N3_consecutive_differenceKL{a});
    Green_rec_N4_consecutive_differencetotalKL(l)=sum(Green_rec_N4_consecutive_differenceKL{a});
    Red_rec_N_consecutive_differencetotalKL(l)=sum(Red_rec_N_consecutive_differenceKL{a});
    Red_rec_N2_consecutive_differencetotalKL(l)=sum(Red_rec_N2_consecutive_differenceKL{a});
    Red_rec_N3_consecutive_differencetotalKL(l)=sum(Red_rec_N3_consecutive_differenceKL{a});
    Red_rec_N4_consecutive_differencetotalKL(l)=sum(Red_rec_N4_consecutive_differenceKL{a});
    Conversion_rec_N_consecutive_differencetotalKL(l)=sum(Conversion_rec_N_consecutive_differenceKL{a});
    Conversion_rec_N2_consecutive_differencetotalKL(l)=sum(Conversion_rec_N2_consecutive_differenceKL{a});
    Conversion_rec_N3_consecutive_differencetotalKL(l)=sum(Conversion_rec_N3_consecutive_differenceKL{a});
    Conversion_rec_N4_consecutive_differencetotalKL(l)=sum(Conversion_rec_N4_consecutive_differenceKL{a});
    
    a=a+1;
end


%%Normalization of the KL divergence values
GreenpretotalKL_norm=(GreenpretotalKL)./max(GreenpretotalKL);
GreenconsecutivetotalKL_norm=(GreenconsecutivetotalKL)./max(GreenconsecutivetotalKL);
RedpretotalKL_norm=(RedpretotalKL)./max(RedpretotalKL);
RedconsecutivetotalKL_norm=(RedconsecutivetotalKL)./max(RedconsecutivetotalKL);
ConversionpretotalKL_norm=(ConversionpretotalKL)./max(ConversionpretotalKL);
ConversionconsecutivetotalKL_norm=(ConversionconsecutivetotalKL)./max(ConversionconsecutivetotalKL);

Green_NpretotalKL_norm=(Green_NpretotalKL)./max(Green_NpretotalKL);
Green_N2pretotalKL_norm=(Green_N2pretotalKL)./max(Green_N2pretotalKL);
Green_N3pretotalKL_norm=(Green_N3pretotalKL)./max(Green_N3pretotalKL);
Green_N4pretotalKL_norm=(Green_N4pretotalKL)./max(Green_N4pretotalKL);
Green_NconsecutivetotalKL_norm=(Green_NconsecutivetotalKL)./max(Green_NconsecutivetotalKL);
Green_N2consecutivetotalKL_norm=(Green_N2consecutivetotalKL)./max(Green_N2consecutivetotalKL);
Green_N3consecutivetotalKL_norm=(Green_N3consecutivetotalKL)./max(Green_N3consecutivetotalKL);
Green_N4consecutivetotalKL_norm=(Green_N4consecutivetotalKL)./max(Green_N4consecutivetotalKL);
Red_NpretotalKL_norm=(Red_NpretotalKL)./max(Red_NpretotalKL);
Red_N2pretotalKL_norm=(Red_N2pretotalKL)./max(Red_N2pretotalKL); 
Red_N3pretotalKL_norm=(Red_N3pretotalKL)./max(Red_N3pretotalKL);
Red_N4pretotalKL_norm=(Red_N4pretotalKL)./max(Red_N4pretotalKL);
Red_NconsecutivetotalKL_norm=(Red_NconsecutivetotalKL)./max(Red_NconsecutivetotalKL); 
Red_N2consecutivetotalKL_norm=(Red_N2consecutivetotalKL)./max(Red_N2consecutivetotalKL);
Red_N3consecutivetotalKL_norm=(Red_N3consecutivetotalKL)./max(Red_N3consecutivetotalKL);
Red_N4consecutivetotalKL_norm=(Red_N4consecutivetotalKL)./max(Red_N4consecutivetotalKL);
Conversion_NpretotalKL_norm=(Conversion_NpretotalKL)./max(Conversion_NpretotalKL);
Conversion_N2pretotalKL_norm=(Conversion_N2pretotalKL)./max(Conversion_N2pretotalKL);
Conversion_N3pretotalKL_norm=(Conversion_N3pretotalKL)./max(Conversion_N3pretotalKL);
Conversion_N4pretotalKL_norm=(Conversion_N4pretotalKL)./max(Conversion_N4pretotalKL);
Conversion_NconsecutivetotalKL_norm=(Conversion_NconsecutivetotalKL)./max(Conversion_NconsecutivetotalKL);
Conversion_N2consecutivetotalKL_norm=(Conversion_N2consecutivetotalKL)./max(Conversion_N2consecutivetotalKL);
Conversion_N3consecutivetotalKL_norm=(Conversion_N3consecutivetotalKL)./max(Conversion_N3consecutivetotalKL);
Conversion_N4consecutivetotalKL_norm=(Conversion_N4consecutivetotalKL)./max(Conversion_N4consecutivetotalKL);

Green_rec_N_pre_differencetotalKL_norm=(Green_rec_N_pre_differencetotalKL)./max(Green_rec_N_pre_differencetotalKL);
Green_rec_N2_pre_differencetotalKL_norm=(Green_rec_N2_pre_differencetotalKL)./max(Green_rec_N2_pre_differencetotalKL);
Green_rec_N3_pre_differencetotalKL_norm=(Green_rec_N3_pre_differencetotalKL)./max(Green_rec_N3_pre_differencetotalKL);
Green_rec_N4_pre_differencetotalKL_norm=(Green_rec_N4_pre_differencetotalKL)./max(Green_rec_N4_pre_differencetotalKL);
Green_rec_N_consecutive_differencetotalKL_norm=(Green_rec_N_consecutive_differencetotalKL)./max(Green_rec_N_consecutive_differencetotalKL);
Green_rec_N2_consecutive_differencetotalKL_norm=(Green_rec_N2_consecutive_differencetotalKL)./max(Green_rec_N2_consecutive_differencetotalKL);
Green_rec_N3_consecutive_differencetotalKL_norm=(Green_rec_N3_consecutive_differencetotalKL)./max(Green_rec_N3_consecutive_differencetotalKL);
Green_rec_N4_consecutive_differencetotalKL_norm=(Green_rec_N4_consecutive_differencetotalKL)./max(Green_rec_N4_consecutive_differencetotalKL);
Red_rec_N_pre_differencetotalKL_norm=(Red_rec_N_pre_differencetotalKL)./max(Red_rec_N_pre_differencetotalKL);
Red_rec_N2_pre_differencetotalKL_norm=(Red_rec_N2_pre_differencetotalKL)./max(Red_rec_N2_pre_differencetotalKL);
Red_rec_N3_pre_differencetotalKL_norm=(Red_rec_N3_pre_differencetotalKL)./max(Red_rec_N3_pre_differencetotalKL);
Red_rec_N4_pre_differencetotalKL_norm=(Red_rec_N4_pre_differencetotalKL)./max(Red_rec_N4_pre_differencetotalKL);
Red_rec_N_consecutive_differencetotalKL_norm=(Red_rec_N_consecutive_differencetotalKL)./max(Red_rec_N_consecutive_differencetotalKL);
Red_rec_N2_consecutive_differencetotalKL_norm=(Red_rec_N2_consecutive_differencetotalKL)./max(Red_rec_N2_consecutive_differencetotalKL);
Red_rec_N3_consecutive_differencetotalKL_norm=(Red_rec_N3_consecutive_differencetotalKL)./max(Red_rec_N3_consecutive_differencetotalKL);
Red_rec_N4_consecutive_differencetotalKL_norm=(Red_rec_N4_consecutive_differencetotalKL)./max(Red_rec_N4_consecutive_differencetotalKL);
Conversion_rec_N_pre_differencetotalKL_norm=(Conversion_rec_N_pre_differencetotalKL)./max(Conversion_rec_N_pre_differencetotalKL);
Conversion_rec_N2_pre_differencetotalKL_norm=(Conversion_rec_N2_pre_differencetotalKL)./max(Conversion_rec_N2_pre_differencetotalKL);
Conversion_rec_N3_pre_differencetotalKL_norm=(Conversion_rec_N3_pre_differencetotalKL)./max(Conversion_rec_N3_pre_differencetotalKL);
Conversion_rec_N4_pre_differencetotalKL_norm=(Conversion_rec_N4_pre_differencetotalKL)./max(Conversion_rec_N4_pre_differencetotalKL);
Conversion_rec_N_consecutive_differencetotalKL_norm=(Conversion_rec_N_consecutive_differencetotalKL)./max(Conversion_rec_N_consecutive_differencetotalKL);
Conversion_rec_N2_consecutive_differencetotalKL_norm=(Conversion_rec_N2_consecutive_differencetotalKL)./max(Conversion_rec_N2_consecutive_differencetotalKL);
Conversion_rec_N3_consecutive_differencetotalKL_norm=(Conversion_rec_N3_consecutive_differencetotalKL)./max(Conversion_rec_N3_consecutive_differencetotalKL);
Conversion_rec_N4_consecutive_differencetotalKL_norm=(Conversion_rec_N4_consecutive_differencetotalKL)./max(Conversion_rec_N4_consecutive_differencetotalKL);

Green_rec_N_correspondingtotalKL_norm=(Green_rec_N_correspondingtotalKL)./max(Green_rec_N_correspondingtotalKL);
Green_rec_N2_correspondingtotalKL_norm=(Green_rec_N2_correspondingtotalKL)./max(Green_rec_N2_correspondingtotalKL);
Green_rec_N3_correspondingtotalKL_norm=(Green_rec_N3_correspondingtotalKL)./max(Green_rec_N3_correspondingtotalKL);
Green_rec_N4_correspondingtotalKL_norm=(Green_rec_N4_correspondingtotalKL)./max(Green_rec_N4_correspondingtotalKL);
Red_rec_N_correspondingtotalKL_norm=(Red_rec_N_correspondingtotalKL)./max(Red_rec_N_correspondingtotalKL);
Red_rec_N2_correspondingtotalKL_norm=(Red_rec_N2_correspondingtotalKL)./max(Red_rec_N2_correspondingtotalKL);
Red_rec_N3_correspondingtotalKL_norm=(Red_rec_N3_correspondingtotalKL)./max(Red_rec_N3_correspondingtotalKL);
Red_rec_N4_correspondingtotalKL_norm=(Red_rec_N4_correspondingtotalKL)./max(Red_rec_N4_correspondingtotalKL);
Conversion_rec_N_correspondingtotalKL_norm=(Conversion_rec_N_correspondingtotalKL)./max(Conversion_rec_N_correspondingtotalKL);
Conversion_rec_N2_correspondingtotalKL_norm=(Conversion_rec_N2_correspondingtotalKL)./max(Conversion_rec_N2_correspondingtotalKL);
Conversion_rec_N3_correspondingtotalKL_norm=(Conversion_rec_N3_correspondingtotalKL)./max(Conversion_rec_N3_correspondingtotalKL);
Conversion_rec_N4_correspondingtotalKL_norm=(Conversion_rec_N4_correspondingtotalKL)./max(Conversion_rec_N4_correspondingtotalKL);

%% Writing all normalized KL values into one cell array

all_KL_values_normalized{1}=GreenpretotalKL_norm;
all_KL_values_normalized{2}=GreenconsecutivetotalKL_norm;
all_KL_values_normalized{3}=RedpretotalKL_norm;
all_KL_values_normalized{4}=RedconsecutivetotalKL_norm;
all_KL_values_normalized{5}=ConversionpretotalKL_norm;
all_KL_values_normalized{6}=ConversionconsecutivetotalKL_norm;

all_KL_values_normalized{7}=Green_NpretotalKL_norm;
all_KL_values_normalized{8}=Green_N2pretotalKL_norm;
all_KL_values_normalized{9}=Green_N3pretotalKL_norm;
all_KL_values_normalized{10}=Green_N4pretotalKL_norm;
all_KL_values_normalized{11}=Green_NconsecutivetotalKL_norm;
all_KL_values_normalized{12}=Green_N2consecutivetotalKL_norm;
all_KL_values_normalized{13}=Green_N3consecutivetotalKL_norm;
all_KL_values_normalized{14}=Green_N4consecutivetotalKL_norm;
all_KL_values_normalized{15}=Red_NpretotalKL_norm;
all_KL_values_normalized{16}=Red_N2pretotalKL_norm; 
all_KL_values_normalized{17}=Red_N3pretotalKL_norm;
all_KL_values_normalized{18}=Red_N4pretotalKL_norm;
all_KL_values_normalized{19}=Red_NconsecutivetotalKL_norm; 
all_KL_values_normalized{20}=Red_N2consecutivetotalKL_norm;
all_KL_values_normalized{21}=Red_N3consecutivetotalKL_norm;
all_KL_values_normalized{22}=Red_N4consecutivetotalKL_norm;
all_KL_values_normalized{23}=Conversion_NpretotalKL_norm;
all_KL_values_normalized{24}=Conversion_N2pretotalKL_norm;
all_KL_values_normalized{25}=Conversion_N3pretotalKL_norm;
all_KL_values_normalized{26}=Conversion_N4pretotalKL_norm;
all_KL_values_normalized{27}=Conversion_NconsecutivetotalKL_norm;
all_KL_values_normalized{28}=Conversion_N2consecutivetotalKL_norm;
all_KL_values_normalized{29}=Conversion_N3consecutivetotalKL_norm;
all_KL_values_normalized{30}=Conversion_N4consecutivetotalKL_norm;

all_KL_values_normalized{31}=Green_rec_N_pre_differencetotalKL_norm;
all_KL_values_normalized{32}=Green_rec_N2_pre_differencetotalKL_norm;
all_KL_values_normalized{33}=Green_rec_N3_pre_differencetotalKL_norm;
all_KL_values_normalized{34}=Green_rec_N4_pre_differencetotalKL_norm;
all_KL_values_normalized{35}=Green_rec_N_consecutive_differencetotalKL_norm;
all_KL_values_normalized{36}=Green_rec_N2_consecutive_differencetotalKL_norm;
all_KL_values_normalized{37}=Green_rec_N3_consecutive_differencetotalKL_norm;
all_KL_values_normalized{38}=Green_rec_N4_consecutive_differencetotalKL_norm;
all_KL_values_normalized{39}=Red_rec_N_pre_differencetotalKL_norm;
all_KL_values_normalized{40}=Red_rec_N2_pre_differencetotalKL_norm;
all_KL_values_normalized{41}=Red_rec_N3_pre_differencetotalKL_norm;
all_KL_values_normalized{42}=Red_rec_N4_pre_differencetotalKL_norm;
all_KL_values_normalized{43}=Red_rec_N_consecutive_differencetotalKL_norm;
all_KL_values_normalized{44}=Red_rec_N2_consecutive_differencetotalKL_norm;
all_KL_values_normalized{45}=Red_rec_N3_consecutive_differencetotalKL_norm;
all_KL_values_normalized{46}=Red_rec_N4_consecutive_differencetotalKL_norm;
all_KL_values_normalized{47}=Conversion_rec_N_pre_differencetotalKL_norm;
all_KL_values_normalized{48}=Conversion_rec_N2_pre_differencetotalKL_norm;
all_KL_values_normalized{49}=Conversion_rec_N3_pre_differencetotalKL_norm;
all_KL_values_normalized{50}=Conversion_rec_N4_pre_differencetotalKL_norm;
all_KL_values_normalized{51}=Conversion_rec_N_consecutive_differencetotalKL_norm;
all_KL_values_normalized{52}=Conversion_rec_N2_consecutive_differencetotalKL_norm;
all_KL_values_normalized{53}=Conversion_rec_N3_consecutive_differencetotalKL_norm;
all_KL_values_normalized{54}=Conversion_rec_N4_consecutive_differencetotalKL_norm;

all_KL_values_normalized{55}=Green_rec_N_correspondingtotalKL_norm;
all_KL_values_normalized{56}=Green_rec_N2_correspondingtotalKL_norm;
all_KL_values_normalized{57}=Green_rec_N3_correspondingtotalKL_norm;
all_KL_values_normalized{58}=Green_rec_N4_correspondingtotalKL_norm;
all_KL_values_normalized{59}=Red_rec_N_correspondingtotalKL_norm;
all_KL_values_normalized{60}=Red_rec_N2_correspondingtotalKL_norm;
all_KL_values_normalized{61}=Red_rec_N3_correspondingtotalKL_norm;
all_KL_values_normalized{62}=Red_rec_N4_correspondingtotalKL_norm;
all_KL_values_normalized{63}=Conversion_rec_N_correspondingtotalKL_norm;
all_KL_values_normalized{64}=Conversion_rec_N2_correspondingtotalKL_norm;
all_KL_values_normalized{65}=Conversion_rec_N3_correspondingtotalKL_norm;
all_KL_values_normalized{66}=Conversion_rec_N4_correspondingtotalKL_norm;

%% writing all raw KL values into one cell array 

all_KL_values{1}=GreenpretotalKL;
all_KL_values{2}=GreenconsecutivetotalKL;
all_KL_values{3}=RedpretotalKL;
all_KL_values{4}=RedconsecutivetotalKL;
all_KL_values{5}=ConversionpretotalKL;
all_KL_values{6}=ConversionconsecutivetotalKL;

all_KL_values{7}=Green_NpretotalKL;
all_KL_values{8}=Green_N2pretotalKL;
all_KL_values{9}=Green_N3pretotalKL;
all_KL_values{10}=Green_N4pretotalKL;
all_KL_values{11}=Green_NconsecutivetotalKL;
all_KL_values{12}=Green_N2consecutivetotalKL;
all_KL_values{13}=Green_N3consecutivetotalKL;
all_KL_values{14}=Green_N4consecutivetotalKL;
all_KL_values{15}=Red_NpretotalKL;
all_KL_values{16}=Red_N2pretotalKL; 
all_KL_values{17}=Red_N3pretotalKL;
all_KL_values{18}=Red_N4pretotalKL;
all_KL_values{19}=Red_NconsecutivetotalKL;
all_KL_values{20}=Red_N2consecutivetotalKL;
all_KL_values{21}=Red_N3consecutivetotalKL;
all_KL_values{22}=Red_N4consecutivetotalKL;
all_KL_values{23}=Conversion_NpretotalKL;
all_KL_values{24}=Conversion_N2pretotalKL;
all_KL_values{25}=Conversion_N3pretotalKL;
all_KL_values{26}=Conversion_N4pretotalKL;
all_KL_values{27}=Conversion_NconsecutivetotalKL;
all_KL_values{28}=Conversion_N2consecutivetotalKL;
all_KL_values{29}=Conversion_N3consecutivetotalKL;
all_KL_values{30}=Conversion_N4consecutivetotalKL;

all_KL_values{31}=Green_rec_N_pre_differencetotalKL;
all_KL_values{32}=Green_rec_N2_pre_differencetotalKL;
all_KL_values{33}=Green_rec_N3_pre_differencetotalKL;
all_KL_values{34}=Green_rec_N4_pre_differencetotalKL;
all_KL_values{35}=Green_rec_N_consecutive_differencetotalKL;
all_KL_values{36}=Green_rec_N2_consecutive_differencetotalKL;
all_KL_values{37}=Green_rec_N3_consecutive_differencetotalKL;
all_KL_values{38}=Green_rec_N4_consecutive_differencetotalKL;
all_KL_values{39}=Red_rec_N_pre_differencetotalKL;
all_KL_values{40}=Red_rec_N2_pre_differencetotalKL;
all_KL_values{41}=Red_rec_N3_pre_differencetotalKL;
all_KL_values{42}=Red_rec_N4_pre_differencetotalKL;
all_KL_values{43}=Red_rec_N_consecutive_differencetotalKL;
all_KL_values{44}=Red_rec_N2_consecutive_differencetotalKL;
all_KL_values{45}=Red_rec_N3_consecutive_differencetotalKL;
all_KL_values{46}=Red_rec_N4_consecutive_differencetotalKL;
all_KL_values{47}=Conversion_rec_N_pre_differencetotalKL;
all_KL_values{48}=Conversion_rec_N2_pre_differencetotalKL;
all_KL_values{49}=Conversion_rec_N3_pre_differencetotalKL;
all_KL_values{50}=Conversion_rec_N4_pre_differencetotalKL;
all_KL_values{51}=Conversion_rec_N_consecutive_differencetotalKL;
all_KL_values{52}=Conversion_rec_N2_consecutive_differencetotalKL;
all_KL_values{53}=Conversion_rec_N3_consecutive_differencetotalKL;
all_KL_values{54}=Conversion_rec_N4_consecutive_differencetotalKL;

all_KL_values{55}=Green_rec_N_correspondingtotalKL;
all_KL_values{56}=Green_rec_N2_correspondingtotalKL;
all_KL_values{57}=Green_rec_N3_correspondingtotalKL;
all_KL_values{58}=Green_rec_N4_correspondingtotalKL;
all_KL_values{59}=Red_rec_N_correspondingtotalKL;
all_KL_values{60}=Red_rec_N2_correspondingtotalKL;
all_KL_values{61}=Red_rec_N3_correspondingtotalKL;
all_KL_values{62}=Red_rec_N4_correspondingtotalKL;
all_KL_values{63}=Conversion_rec_N_correspondingtotalKL;
all_KL_values{64}=Conversion_rec_N2_correspondingtotalKL;
all_KL_values{65}=Conversion_rec_N3_correspondingtotalKL;
all_KL_values{66}=Conversion_rec_N4_correspondingtotalKL;


%% Curve fitting to the raw KL values
%setting up fit type and options, smoothing spline is chosen in order to
%not lose too many information, smoothing parameter is chosen relatively
%close to 1 which results in rougher curve but more accurate data info
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.95;

%fitting model separately to each totalKL found
s=[10 20 50 100 250 450 650 850];
sc=[0 10 20 50 100 250 450 650 850];
scspikes=(0:1:850);
scspikes_vertical=vertcat(scspikes);
spikes=(10:1:850);
spikes_vertical=vertcat(spikes);
%prepares data to fit the curve, arranges data horizontally etc.
for i=1:size(all_KL_values,2)
    if size(all_KL_values{i},2)==9
        j=sc; 
        spikes_vertical_data=scspikes_vertical;
    else
        j=s;spikes_vertical_data=spikes_vertical;
    end
   %y=horzcat(all_KL_values{i})    
[xData, yData] = prepareCurveData( j, all_KL_values{i});
[fitresult{i}, gof{i}] = fit( xData, yData, ft, opts);
[fx{i}, fxx{i}]=differentiate(fitresult{i},spikes_vertical_data);
fitresult_fitted_y_values{i}=feval(fitresult{i},xData);
fastest_change_y_value{i}=max(fx{i});
idx_max=(fx{i}==fastest_change_y_value{i});
fastest_change_spike_value{i}=spikes_vertical_data(idx_max);

%%determining the dynamic range 
% the largest range in which the first derivative stays postive (the widest range that the original function is increasing)
% there are experiments where the conversion starts early and even though
% as a result of fit the dfunction seems decreasing during intervals, but
% this could be disgarded because small spike conversions are not
% consistent througout all experiments.
subplot(3,1,1)
plot( fitresult{i}, xData, yData );
subplot(3,1,2)
plot(spikes_vertical_data,fx{i});
subplot(3,1,3)
plot(spikes_vertical_data,fxx{i});

sign_fx=sign(fx{i});
idx_positive=(sign_fx==1);
%extracting the biggest uninterrupted positive region
o=1;
z=1;
    for k=1:size(idx_positive,1)
        if k<size(idx_positive,1) & idx_positive((k+1),1)~=idx_positive(k,1) & idx_positive((k+1),1)==1
            change_from_zero_to_one{i}(o)=k;
            o=o+1;
     elseif k<size(idx_positive,1) & idx_positive((k+1),1)~=idx_positive(k,1) & idx_positive((k+1),1)==0
            change_from_one_to_zero{i}(z)=k;
            z=z+1;
        end
    end
% change_from_zero_to_one=horzcat(change_from_zero_to_one{:});
% change_from_one_to_zero=horzcat(change_from_one_to_zero{:});
% last_zero_to_one_change=max(change_from_zero_to_one{i});
% last_one_to_zero_change=max(change_from_one_to_zero{i});
% if last_zero_to_one_change>last_one_to_zero_change
%     dynamic_range{i}=[(last_zero_to_one_change+10),850];
% else
%     dynamic_range{i}=[(last_zero_to_one_change),last_one_to_zero_change];
% end
clearvars last_zero_to_one_change last_one_to_zero_change idx_positive sign_fx idx_max y xData yData
end




t=1








% % Set up fittype and options.
% ft = fittype( 'pchipinterp' );
% s=[10 20 50 100 250 450 650 850];
% 
% [x, yData{1}] = prepareCurveData(s, all_KL_values{1});
% [fitresult{1}, gof{1}] = fit(x, yData{1}, ft , 'Normalize', 'on');
% a=zeros(1,7);
% a=coeffvalues(fitresult{1});
% syms k;
% fitresult_formula{1}=poly2sym(a,k);
% 
% fitresult_firstderivative{1}=diff(fitresult_formula{1});
% fitresult_firstderivative_all_y_values{1}=subs(fitresult_firstderivative{1},(10:1:850));
% [q, zData{1}] = prepareCurveData((10:1:850),fitresult_firstderivative_all_y_values{1});
% [fitresult{2}, gof{1}]=fit(q, zData{1}, ft, 'Normalize', 'on');
% 
% plot(fitresult{2}, q, zData{1});
% 
% fitresult_firstderivative_fastest_change_y_value{1}=max(fitresult_firstderivative_all_y_values{1});
% % find the corresponding tspike values to the y value
% assume(k, 'real');
% assumeAlso(k>0);
% fitresult_firstderivative_fastest_change_spike_value{1}=solve(fitresult_firstderivative{1}==fitresult_firstderivative_fastest_change_y_value{1},k);
% 
% 
% fitresult_secondderivative{1}=diff(fitresult_formula{1},2);
% fitresult_secondderivative_all_y_values{1}=subs(fitresult_secondderivative{1},(10:1:850));
% [q, zData{1}] = prepareCurveData((10:1:850),fitresult_secondderivative_all_y_values{1});
% [fitresult{2}, gof{1}]=fit(q, zData{1}, ft, 'Normalize', 'on');
% plot(fitresult{2}, q, zData{1});
% 
% 
% %%second derivative part
% fitresult_secondderivative{1}=diff(fitresult_formula{1},2); % find the maximum and minimum s values where second derivative is equal to zero, this will give us the dynamic range for that cell in which KL is effected by the spikes.
% assume(k,'real');
% assumeAlso(k>0);
% fitresult_secondderivative_all_spike_values{1}=solve(fitresult_secondderivative{1}==0,k);
% %min of these values is the spike count where conversion starts, and the
% %max value is the spike value where conversion ends
% fitresult_secondderivative_all_spike_values{1}=double(fitresult_secondderivative_all_spike_values{1});
% [q, zData{1}] = prepareCurveData((10:1:850),fitresult_secondderivative_all_spike_values{1});
% [fitresult{2}, gof{1}]=fit(q, zData{1}, ft , 'Normalize', 'on');
% 
% plot(fitresult{2}, q, zData{1});
% fitresult_secondderivative_change_start_spike_value{1}=min(fitresult_secondderivative_all_spike_values{1});
% fitresult_secondderivative_change_end_spike_value{1}=max(fitresult_secondderivative_all_spike_values{1});
% 
% 
% for i=1:66
%     syms x;
% [x, yData{i}] = prepareCurveData(s, all_KL_values{i});
% % Fit model to data.
% 
% [fitresult{i}, gof{i}] = fit(x, yData{i}, ft, 'Normalize', 'on' );
% 
% 
% fitresult_firstderivative{i}=diff(fitresult{i}); %find the s value where first derivative is highest to determine the fastest change point
% fitresult_firstderivative_all_y_values{i}=solve(fitresult_firstderivative{i},(10:1:850));
% fitresult_firstderivative_fastest_change_y_value{i}=max(fitresult_firstderivative_all_y_values{i});
% % find the corresponding tspike values to the y value
% fitresult_firstderivative_fastest_change_spike_value{i}=solve(fitresult_firstderivative{i},fitresult_firstderivative_fastest_change_y_value{i});
% 
% fitresult_secondderivative{i}=diff(fitresult{i},2); % find the maximum and minimum s values where second derivative is equal to zero, this will give us the dynamic range for that cell in which KL is effected by the spikes.
% fitresult_secondderivative_all_spike_values{i}=solve(fitresult_secondderivative{i}==0,x,AllRealRoots);
% %min of these values is the spike count where conversion starts, and the
% %max value is the spike value where conversion ends
% fitresult_secondderivative_change_start_spike_value{i}=min(fitresult_secondderivative_all_spike_values{i});
% fitresult_secondderivative_change_end_spike_value{i}=max(fitresult_secondderivative_all_spike_values{i});
% 
% end

% figure('Name','Raw KL data with fitted curve');
% plot( fitresult{i}, s, yData{i} );
% % Generate 5 hue-saturation-value color map for your data
% colorVec = ['g-*','g-o','k-*','k-o','m-*','m-o',
% 'r-*','r-o','k-*','k-o','m-*','m-o',
% 'b-*','b-o','k-*','k-o','m-*','m-o'];
% % Plot and change the color for each line
% hold on;
% for i = 1:fitresult{i}
%     subplot(5,3,i);
%     plot(fitresult{i}, s, yData{i},'Color',colorVec(i));
% end
% hold off;
% 
% hold on;
% for i = 1:fitresult{i}
%     subplot(5,3,i);
%     plot(fitresult{i}, s, yData{i},'Color',colorVec(i));
% end
% hold off;
% 
% 
% % Set up fittype and options.
% ft = fittype( 'poly6' ); 
% 
% % Fit model to data.
% 
% [fitresult, gof] = fit( s, yData, ft );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, s, yData );
% legend( h, 'M_red_consecutive vs. xmatrix', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel xmatrix
% ylabel M_red_consecutive
% grid on

%% Curve fitting to the normalized KL values

% Set up fittype and options.
% 
% for i=1:66
%     syms x;
% [x, yData_normalized{i}] = prepareCurveData(s, all_KL_values_normalized{i});
% % Fit model to data.
% 
% [fitresult_normalized{i}, gof_normalized{i}] = fit(x, yData_normalized{i}, ft );
% 
% fitresult_firstderivative_normalized{i}=diff(fitresult_normalized{i}); %find the s value where first derivative is highest to determine the fastest change point
% fitresult_firstderivative_all_y_values_normalized{i}=solve(fitresult_firstderivative_normalized{i},(10:1:850));
% fitresult_firstderivative_fastest_change_y_value_normalized{i}=max(fitresult_firstderivative_all_y_values_normalized{i});
% % find the corresponding tspike values to the y value
% fitresult_firstderivative_fastest_change_spike_value_normalized{i}=solve(fitresult_firstderivative_normalized{i},fitresult_firstderivative_fastest_change_y_value_normalized{i});
% 
% fitresult_secondderivative_normalized{i}=diff(fitresult_normalized{i},2); % find the maximum and minimum s values where second derivative is equal to zero, this will give us the dynamic range for that cell in which KL is effected by the spikes.
% fitresult_secondderivative_allvalues_normalized{i}=solve(fitresult_secondderivative_normalized{i}==0,x,AllRealRoots);
% %min of these values is the spike count where conversion starts, and the
% %max value is the spike value where conversion ends
% fitresult_secondderivative_change_start_spike_value_normalized{i}=min(fitresult_secondderivative_allvalues_normalized{i});
% fitresult_secondderivative_change_end_spike_value_normalized{i}=max(fitresult_secondderivative_allvalues_normalized{i});
% 
% end

% figure('Name','Normalized KL data with fitted curve');
% plot( all_KL_values_normalized{i}, s, yData_normalized{i} );
% % Generate 5 hue-saturation-value color map for your data
% colorVec = ['g-*','g-o','k-*','k-o','m-*','m-o',
% 'r-*','r-o','k-*','k-o','m-*','m-o',
% 'b-*','b-o','k-*','k-o','m-*','m-o'];
% % Plot and change the color for each line
% hold on;
% for i = 1:fitresult{i}
%     subplot(5,3,i);
%     plot(fitresult{i}, s, yData{i},'Color',colorVec(i));
% end
% hold off;
% 
% hold on;
% for i = 1:fitresult{i}
%     subplot(5,3,i);
%     plot(fitresult{i}, s, yData{i},'Color',colorVec(i));
% end
% hold off;
% 
% % Set up fittype and options.
% ft = fittype( 'poly6' ); 
% 
% % Fit model to data.
% 
% [fitresult, gof] = fit( s, yData, ft );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, s, yData );
% legend( h, 'M_red_consecutive vs. xmatrix', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel xmatrix
% ylabel M_red_consecutive
% grid on
% 




       
%% Writing the raw calculated KL values to a file in order to analyze all experiments together
if ~(nargin<1)

% writing the recorded neuron KL divergence values
filename1=('Recorded neurons KL divergence values 850spikes.csv');
M=csvread(filename1);
Mfinal=[M; GreenpretotalKL; GreenconsecutivetotalKL; 
RedpretotalKL; RedconsecutivetotalKL; 
ConversionpretotalKL; ConversionconsecutivetotalKL];

csvwrite(filename1,Mfinal);
Mfinal;

% wiriting the non recorded neurons KL divergence values

filename2=('Non-recorded neurons KL divergence values 850spikes.csv');
S=csvread(filename2);
Sfinal=[S; Green_NpretotalKL; Green_N2pretotalKL; Green_N3pretotalKL;Green_N4pretotalKL;
Green_NconsecutivetotalKL; Green_N2consecutivetotalKL; Green_N3consecutivetotalKL; Green_N4consecutivetotalKL; 
Red_NpretotalKL; Red_N2pretotalKL; Red_N3pretotalKL; Red_N4pretotalKL; 
Red_NconsecutivetotalKL; Red_N2consecutivetotalKL; Red_N3consecutivetotalKL; Red_N4consecutivetotalKL;
Conversion_NpretotalKL; Conversion_N2pretotalKL; Conversion_N3pretotalKL; Conversion_N4pretotalKL; 
Conversion_NconsecutivetotalKL; Conversion_N2consecutivetotalKL; Conversion_N3consecutivetotalKL; 
Conversion_N4consecutivetotalKL];

csvwrite(filename2,Sfinal);
Sfinal;

% wiriting the non recorded-recorded neurons KL divergence values

filename3=('Non-recorded neurons compared recorded neurons KL divergence values 850spikes.csv');
T=csvread(filename3);
Tfinal=[T;Green_rec_N_pre_differencetotalKL; Green_rec_N2_pre_differencetotalKL; Green_rec_N3_pre_differencetotalKL;
Green_rec_N4_pre_differencetotalKL; 
Green_rec_N_consecutive_differencetotalKL; Green_rec_N2_consecutive_differencetotalKL;Green_rec_N3_consecutive_differencetotalKL; 
Green_rec_N4_consecutive_differencetotalKL;
Red_rec_N_pre_differencetotalKL;Red_rec_N2_pre_differencetotalKL; Red_rec_N3_pre_differencetotalKL;
Red_rec_N4_pre_differencetotalKL; 
Red_rec_N_consecutive_differencetotalKL; Red_rec_N2_consecutive_differencetotalKL; Red_rec_N3_consecutive_differencetotalKL
Red_rec_N4_consecutive_differencetotalKL; 
Conversion_rec_N_pre_differencetotalKL; Conversion_rec_N2_pre_differencetotalKL; Conversion_rec_N3_pre_differencetotalKL; 
Conversion_rec_N4_pre_differencetotalKL;
Conversion_rec_N_consecutive_differencetotalKL;Conversion_rec_N2_consecutive_differencetotalKL; Conversion_rec_N3_consecutive_differencetotalKL;
Conversion_rec_N4_consecutive_differencetotalKL];

csvwrite(filename3,Tfinal);
Tfinal;

% wiriting the non recorded-recorded neurons corresponding condition KL
% divergence values, it is written to separate file from the other recorded
% nonrecorded comparisons because the lenght of the vector for
% corresponding condition is 10 instead of 9 like the others


filename4=('Non-recorded neurons compared recorded neurons corresponding condition KL divergence values 850spikes.csv');
Z=csvread(filename4);
Zfinal=[Z; Green_rec_N_correspondingtotalKL; Green_rec_N2_correspondingtotalKL; Green_rec_N3_correspondingtotalKL;
Green_rec_N4_correspondingtotalKL;
Red_rec_N_correspondingtotalKL; Red_rec_N2_correspondingtotalKL; Red_rec_N3_correspondingtotalKL;
Red_rec_N4_correspondingtotalKL;
Conversion_rec_N_correspondingtotalKL; Conversion_rec_N2_correspondingtotalKL; Conversion_rec_N3_correspondingtotalKL;
Conversion_rec_N4_correspondingtotalKL];


csvwrite(filename4,Zfinal);
Zfinal;

%writing the distances 
filename5=('Nonrecorded neuron distances 850spikes.csv');
P=csvread(filename5);
Pfinal=[P; nonrecorded_distance; nonrecorded2_distance; nonrecorded3_distance;
nonrecorded4_distance];

csvwrite(filename5,Pfinal);

%% Writing the normalized calculated KL values to a file in order to analyze all experiments together


% writing the recorded neuron KL divergence values
filename6=('Normalized Recorded neurons KL divergence values 850spikes.csv');
A=csvread(filename6);
Afinal=[A; GreenpretotalKL_norm; GreenconsecutivetotalKL_norm; 
RedpretotalKL_norm; RedconsecutivetotalKL_norm; 
ConversionpretotalKL_norm; ConversionconsecutivetotalKL_norm];

csvwrite(filename6,Afinal);
Afinal;

% wiriting the normalized non recorded neurons KL divergence values

filename7=('Normalized Non-recorded neurons KL divergence values 850spikes.csv');
B=csvread(filename7);
Bfinal=[B; Green_NpretotalKL_norm; Green_N2pretotalKL_norm; Green_N3pretotalKL_norm;Green_N4pretotalKL_norm;
Green_NconsecutivetotalKL_norm; Green_N2consecutivetotalKL_norm; Green_N3consecutivetotalKL_norm; Green_N4consecutivetotalKL_norm; 
Red_NpretotalKL_norm; Red_N2pretotalKL; Red_N3pretotalKL_norm; Red_N4pretotalKL_norm; 
Red_NconsecutivetotalKL_norm; Red_N2consecutivetotalKL_norm; Red_N3consecutivetotalKL_norm; Red_N4consecutivetotalKL_norm;
Conversion_NpretotalKL_norm; Conversion_N2pretotalKL_norm; Conversion_N3pretotalKL_norm; Conversion_N4pretotalKL_norm; 
Conversion_NconsecutivetotalKL_norm; Conversion_N2consecutivetotalKL_norm; Conversion_N3consecutivetotalKL_norm; 
Conversion_N4consecutivetotalKL_norm];

csvwrite(filename7,Bfinal);
Bfinal;

% wiriting the normalized non recorded-recorded neurons KL divergence values

filename8=('Normalized Non-recorded neurons compared recorded neurons KL divergence values 850spikes.csv');
C=csvread(filename8);
Cfinal=[C;Green_rec_N_pre_differencetotalKL_norm; Green_rec_N2_pre_differencetotalKL_norm; Green_rec_N3_pre_differencetotalKL_norm;
Green_rec_N4_pre_differencetotalKL_norm; 
Green_rec_N_consecutive_differencetotalKL_norm; Green_rec_N2_consecutive_differencetotalKL_norm;Green_rec_N3_consecutive_differencetotalKL_norm; 
Green_rec_N4_consecutive_differencetotalKL_norm;
Red_rec_N_pre_differencetotalKL_norm;Red_rec_N2_pre_differencetotalKL_norm; Red_rec_N3_pre_differencetotalKL_norm;
Red_rec_N4_pre_differencetotalKL_norm; 
Red_rec_N_consecutive_differencetotalKL_norm; Red_rec_N2_consecutive_differencetotalKL_norm; Red_rec_N3_consecutive_differencetotalKL_norm;
Red_rec_N4_consecutive_differencetotalKL_norm; 
Conversion_rec_N_pre_differencetotalKL_norm; Conversion_rec_N2_pre_differencetotalKL_norm; Conversion_rec_N3_pre_differencetotalKL_norm; 
Conversion_rec_N4_pre_differencetotalKL_norm;
Conversion_rec_N_consecutive_differencetotalKL_norm;Conversion_rec_N2_consecutive_differencetotalKL_norm; Conversion_rec_N3_consecutive_differencetotalKL_norm;
Conversion_rec_N4_consecutive_differencetotalKL_norm];

csvwrite(filename8,Cfinal);
Cfinal;

% wiriting the normalized non recorded-recorded neurons corresponding condition KL
% divergence values, it is written to separate file from the other recorded
% nonrecorded comparisons because the length of the vector for
% corresponding condition is 10 instead of 9 like the others


filename9=('Normalized Non-recorded neurons compared recorded neurons corresponding condition KL divergence values 850spikes.csv');
E=csvread(filename9);
Efinal=[E; Green_rec_N_correspondingtotalKL_norm; Green_rec_N2_correspondingtotalKL_norm; Green_rec_N3_correspondingtotalKL_norm;
Green_rec_N4_correspondingtotalKL_norm;
Red_rec_N_correspondingtotalKL_norm; Red_rec_N2_correspondingtotalKL_norm; Red_rec_N3_correspondingtotalKL_norm;
Red_rec_N4_correspondingtotalKL_norm;
Conversion_rec_N_correspondingtotalKL_norm; Conversion_rec_N2_correspondingtotalKL_norm; Conversion_rec_N3_correspondingtotalKL_norm;
Conversion_rec_N4_correspondingtotalKL_norm];


csvwrite(filename9,Efinal);
Efinal;

%% writing fitresult values into csv files
%% writing fitresult_normalized values into csv files



%% Plotting all the raw KL divergences

%spike conditions
s=[10 20 50 100 250 450 650 850];
sc=[0 10 20 50 100 250 450 650 850];

%% First non recorded
%Green recorded and first non recorded neuron plots
figure('Name','Raw KL data without fitted curve');
subplot(5,3,1);
plot(s,GreenpretotalKL,'g-*',s,GreenconsecutivetotalKL,'g-o',s,Green_NpretotalKL,'k-*',s,Green_NconsecutivetotalKL,'k-o',s,Green_rec_N_pre_differencetotalKL,'m-*',s,Green_rec_N_consecutive_differencetotalKL,'m-o');
title('Green channel, recorded, first nonrecorded, recorded-first nonrecorded difference KL divergence values');
%Red recorded and first non recorded neuron plots
subplot(5,3,2);
plot(s,RedpretotalKL,'r-*',s,RedconsecutivetotalKL,'r-o',s,Red_NpretotalKL,'k-*',s,Red_NconsecutivetotalKL,'k-o',s,Red_rec_N_pre_differencetotalKL,'m-*',s,Red_rec_N_consecutive_differencetotalKL,'m-o');
title('Red channel, recorded, first nonrecorded, recorded-first nonrecorded difference KL divergence values');
%Conversion rate recorded and first non recorded neuron plots
subplot(5,3,3);
plot(s,ConversionpretotalKL,'b-*',s,ConversionconsecutivetotalKL,'b-o',s,Conversion_NpretotalKL,'k-*',s,Conversion_NconsecutivetotalKL,'k-o',s,Conversion_rec_N_pre_differencetotalKL,'m-*',s,Conversion_rec_N_consecutive_differencetotalKL,'m-o');
title('Conversion rate, recorded, first nonrecorded, recorded-first nonrecorded difference KL divergence values');

%% Second non recorded
%Green recorded and second non recorded neuron plots
subplot(5,3,4);
plot(s,GreenpretotalKL,'g-*',s,GreenconsecutivetotalKL,'g-o',s,Green_N2pretotalKL,'k-*',s,Green_N2consecutivetotalKL,'k-o',s,Green_rec_N2_pre_differencetotalKL,'m-*',s,Green_rec_N2_consecutive_differencetotalKL,'m-o');
title('Green channel, recorded, second nonrecorded, recorded-second nonrecorded difference KL divergence values');
%Red recorded and second non recorded neuron plots
subplot(5,3,5);
plot(s,RedpretotalKL,'r-*',s,RedconsecutivetotalKL,'r-o',s,Red_N2pretotalKL,'k-*',s,Red_N2consecutivetotalKL,'k-o',s,Red_rec_N2_pre_differencetotalKL,'m-*',s,Red_rec_N2_consecutive_differencetotalKL,'m-o');
title('Red channel, recorded, second nonrecorded, recorded-second nonrecorded difference KL divergence values');
%Conversion rate recorded and second non recorded neuron plots
subplot(5,3,6);
plot(s,ConversionpretotalKL,'b-*',s,ConversionconsecutivetotalKL,'b-o',s,Conversion_N2pretotalKL,'k-*',s,Conversion_N2consecutivetotalKL,'k-o',s,Conversion_rec_N2_pre_differencetotalKL,'m-*',s,Conversion_rec_N2_consecutive_differencetotalKL,'m-o');
title('Conversion rate, recorded, second nonrecorded, recorded-second nonrecorded difference KL divergence values');

%% Third non recorded
%Green recorded and third non recorded neuron plots
subplot(5,3,7);
plot(s,GreenpretotalKL,'g-*',s,GreenconsecutivetotalKL,'g-o',s,Green_N3pretotalKL,'k-*',s,Green_N3consecutivetotalKL,'k-o',s,Green_rec_N3_pre_differencetotalKL,'m-*',s,Green_rec_N3_consecutive_differencetotalKL,'m-o');
title('Green channel, recorded, third nonrecorded, recorded-third nonrecorded difference KL divergence values');
%Red recorded and third non recorded neuron plots
subplot(5,3,8);
plot(s,RedpretotalKL,'r-*',s,RedconsecutivetotalKL,'r-o',s,Red_N3pretotalKL,'k-*',s,Red_N3consecutivetotalKL,'k-o',s,Red_rec_N3_pre_differencetotalKL,'m-*',s,Red_rec_N3_consecutive_differencetotalKL,'m-o');
title('Red channel, recorded, third nonrecorded, recorded-third nonrecorded difference KL divergence values');
%Conversion rate recorded and third non recorded neuron plots
subplot(5,3,9);
plot(s,ConversionpretotalKL,'b-*',s,ConversionconsecutivetotalKL,'b-o',s,Conversion_N3pretotalKL,'k-*',s,Conversion_N3consecutivetotalKL,'k-o',s,Conversion_rec_N3_pre_differencetotalKL,'m-*',s,Conversion_rec_N3_consecutive_differencetotalKL,'m-o');
title('Conversion rate, recorded, third nonrecorded, recorded-third nonrecorded difference KL divergence values');
%% Fourth non recorded
%Green recorded and fourth non recorded neuron plots
subplot(5,3,10);
plot(s,GreenpretotalKL,'g-*',s,GreenconsecutivetotalKL,'g-o',s,Green_N4pretotalKL,'k-*',s,Green_N4consecutivetotalKL,'k-o',s,Green_rec_N4_pre_differencetotalKL,'m-*',s,Green_rec_N4_consecutive_differencetotalKL,'m-o');
title('Green channel, recorded, fourth nonrecorded, recorded-fourth nonrecorded difference KL divergence values');
%Red recorded and fourth non recorded neuron plots
subplot(5,3,11);
plot(s,RedpretotalKL,'r-*',s,RedconsecutivetotalKL,'r-o',s,Red_N4pretotalKL,'k-*',s,Red_N4consecutivetotalKL,'k-o',s,Red_rec_N4_pre_differencetotalKL,'m-*',s,Red_rec_N4_consecutive_differencetotalKL,'m-o');
title('Red channel, recorded, fourth nonrecorded, recorded-fourth nonrecorded difference KL divergence values');
%Conversion rate recorded and fourth non recorded neuron plots
subplot(5,3,12);
plot(s,ConversionpretotalKL,'b-*',s,ConversionconsecutivetotalKL,'b-o',s,Conversion_N4pretotalKL,'k-*',s,Conversion_N4consecutivetotalKL,'k-o',s,Conversion_rec_N4_pre_differencetotalKL,'m-*',s,Conversion_rec_N4_consecutive_differencetotalKL,'m-o');
title('Conversion rate, recorded, fourth nonrecorded, recorded-fourth nonrecorded difference KL divergence values');
%% Corresponding condition plotted for each nonrecorded neuron in a separate figure because the lenght is higher then the other comparison conditions
%Green 4 non recorded neurons plot
subplot(5,3,13);
plot(sc,Green_rec_N_correspondingtotalKL,'-ys',sc,Green_rec_N2_correspondingtotalKL,'-ms',sc,Green_rec_N3_correspondingtotalKL,'-cs',sc,Green_rec_N4_correspondingtotalKL,'-ks');
title('Green channel, corresponding condition of all non recorded neurons');
%Red 4 non recorded neurons plot
subplot(5,3,14);
plot(sc,Red_rec_N_correspondingtotalKL,'-ys',sc,Red_rec_N2_correspondingtotalKL,'-ms',sc,Red_rec_N3_correspondingtotalKL,'-cs',sc,Red_rec_N4_correspondingtotalKL,'-ks');
title('Red channel, corresponding condition of all non recorded neurons');
%Conversion rate 4 non recorded neurons plot
subplot(5,3,15);
plot(sc,Conversion_rec_N_correspondingtotalKL,'-ys',sc,Conversion_rec_N2_correspondingtotalKL,'-ms',sc,Conversion_rec_N3_correspondingtotalKL,'-cs',sc,Conversion_rec_N4_correspondingtotalKL,'-ks');
title('Conversion rate, corresponding condition of all non recorded neurons');



save (['KL divergence' finalfile '.mat']);

%% Plotting all the normalized KL divergences

%spike conditions
s=[10 20 50 100 250 450 650 850];
sc=[0 10 20 50 100 250 450 650 850];

%% First non recorded
%Green recorded and first non recorded neuron plots
figure('Name','Normalized KL data without fitted curve');
subplot(5,3,1);
plot(s,GreenpretotalKL_norm,'g-*',s,GreenconsecutivetotalKL_norm,'g-o',s,Green_NpretotalKL_norm,'k-*',s,Green_NconsecutivetotalKL_norm,'k-o',s,Green_rec_N_pre_differencetotalKL_norm,'m-*',s,Green_rec_N_consecutive_differencetotalKL_norm,'m-o');
title('Green channel, recorded, first nonrecorded, recorded-first nonrecorded difference KL divergence values');
%Red recorded and first non recorded neuron plots
subplot(5,3,2);
plot(s,RedpretotalKL_norm,'r-*',s,RedconsecutivetotalKL_norm,'r-o',s,Red_NpretotalKL_norm,'k-*',s,Red_NconsecutivetotalKL_norm,'k-o',s,Red_rec_N_pre_differencetotalKL_norm,'m-*',s,Red_rec_N_consecutive_differencetotalKL_norm,'m-o');
title('Red channel, recorded, first nonrecorded, recorded-first nonrecorded difference KL divergence values');
%Conversion rate recorded and first non recorded neuron plots
subplot(5,3,3);
plot(s,ConversionpretotalKL_norm,'b-*',s,ConversionconsecutivetotalKL_norm,'b-o',s,Conversion_NpretotalKL_norm,'k-*',s,Conversion_NconsecutivetotalKL_norm,'k-o',s,Conversion_rec_N_pre_differencetotalKL_norm,'m-*',s,Conversion_rec_N_consecutive_differencetotalKL_norm,'m-o');
title('Conversion rate, recorded, first nonrecorded, recorded-first nonrecorded difference KL divergence values');

%% Second non recorded
%Green recorded and second non recorded neuron plots
subplot(5,3,4);
plot(s,GreenpretotalKL_norm,'g-*',s,GreenconsecutivetotalKL_norm,'g-o',s,Green_N2pretotalKL_norm,'k-*',s,Green_N2consecutivetotalKL_norm,'k-o',s,Green_rec_N2_pre_differencetotalKL_norm,'m-*',s,Green_rec_N2_consecutive_differencetotalKL_norm,'m-o');
title('Green channel, recorded, second nonrecorded, recorded-second nonrecorded difference KL divergence values');
%Red recorded and second non recorded neuron plots
subplot(5,3,5);
plot(s,RedpretotalKL_norm,'r-*',s,RedconsecutivetotalKL_norm,'r-o',s,Red_N2pretotalKL_norm,'k-*',s,Red_N2consecutivetotalKL_norm,'k-o',s,Red_rec_N2_pre_differencetotalKL_norm,'m-*',s,Red_rec_N2_consecutive_differencetotalKL_norm,'m-o');
title('Red channel, recorded, second nonrecorded, recorded-second nonrecorded difference KL divergence values');
%Conversion rate recorded and second non recorded neuron plots
subplot(5,3,6);
plot(s,ConversionpretotalKL_norm,'b-*',s,ConversionconsecutivetotalKL_norm,'b-o',s,Conversion_N2pretotalKL_norm,'k-*',s,Conversion_N2consecutivetotalKL_norm,'k-o',s,Conversion_rec_N2_pre_differencetotalKL_norm,'m-*',s,Conversion_rec_N2_consecutive_differencetotalKL_norm,'m-o');
title('Conversion rate, recorded, second nonrecorded, recorded-second nonrecorded difference KL divergence values');

%% Third non recorded
%Green recorded and third non recorded neuron plots
subplot(5,3,7);
plot(s,GreenpretotalKL_norm,'g-*',s,GreenconsecutivetotalKL_norm,'g-o',s,Green_N3pretotalKL_norm,'k-*',s,Green_N3consecutivetotalKL_norm,'k-o',s,Green_rec_N3_pre_differencetotalKL_norm,'m-*',s,Green_rec_N3_consecutive_differencetotalKL_norm,'m-o');
title('Green channel, recorded, third nonrecorded, recorded-third nonrecorded difference KL divergence values');
%Red recorded and third non recorded neuron plots
subplot(5,3,8);
plot(s,RedpretotalKL_norm,'r-*',s,RedconsecutivetotalKL_norm,'r-o',s,Red_N3pretotalKL_norm,'k-*',s,Red_N3consecutivetotalKL_norm,'k-o',s,Red_rec_N3_pre_differencetotalKL_norm,'m-*',s,Red_rec_N3_consecutive_differencetotalKL_norm,'m-o');
title('Red channel, recorded, third nonrecorded, recorded-third nonrecorded difference KL divergence values');
%Conversion rate recorded and third non recorded neuron plots
subplot(5,3,9);
plot(s,ConversionpretotalKL_norm,'b-*',s,ConversionconsecutivetotalKL_norm,'b-o',s,Conversion_N3pretotalKL_norm,'k-*',s,Conversion_N3consecutivetotalKL_norm,'k-o',s,Conversion_rec_N3_pre_differencetotalKL_norm,'m-*',s,Conversion_rec_N3_consecutive_differencetotalKL_norm,'m-o');
title('Conversion rate, recorded, third nonrecorded, recorded-third nonrecorded difference KL divergence values');
%% Fourth non recorded
%Green recorded and fourth non recorded neuron plots
subplot(5,3,10);
plot(s,GreenpretotalKL_norm,'g-*',s,GreenconsecutivetotalKL_norm,'g-o',s,Green_N4pretotalKL_norm,'k-*',s,Green_N4consecutivetotalKL_norm,'k-o',s,Green_rec_N4_pre_differencetotalKL_norm,'m-*',s,Green_rec_N4_consecutive_differencetotalKL_norm,'m-o');
title('Green channel, recorded, fourth nonrecorded, recorded-fourth nonrecorded difference KL divergence values');
%Red recorded and fourth non recorded neuron plots
subplot(5,3,11);
plot(s,RedpretotalKL_norm,'r-*',s,RedconsecutivetotalKL_norm,'r-o',s,Red_N4pretotalKL_norm,'k-*',s,Red_N4consecutivetotalKL_norm,'k-o',s,Red_rec_N4_pre_differencetotalKL_norm,'m-*',s,Red_rec_N4_consecutive_differencetotalKL_norm,'m-o');
title('Red channel, recorded, fourth nonrecorded, recorded-fourth nonrecorded difference KL divergence values');
%Conversion rate recorded and fourth non recorded neuron plots
subplot(5,3,12);
plot(s,ConversionpretotalKL_norm,'b-*',s,ConversionconsecutivetotalKL_norm,'b-o',s,Conversion_N4pretotalKL_norm,'k-*',s,Conversion_N4consecutivetotalKL_norm,'k-o',s,Conversion_rec_N4_pre_differencetotalKL_norm,'m-*',s,Conversion_rec_N4_consecutive_differencetotalKL_norm,'m-o');
title('Conversion rate, recorded, fourth nonrecorded, recorded-fourth nonrecorded difference KL divergence values');
%% Corresponding condition plotted for each nonrecorded neuron in a separate figure because the lenght is higher then the other comparison conditions
%Green 4 non recorded neurons plot
subplot(5,3,13);
plot(sc,Green_rec_N_correspondingtotalKL_norm,'-ys',sc,Green_rec_N2_correspondingtotalKL_norm,'-ms',sc,Green_rec_N3_correspondingtotalKL_norm,'-cs',sc,Green_rec_N4_correspondingtotalKL_norm,'-ks');
title('Green channel, corresponding condition of all non recorded neurons');
%Red 4 non recorded neurons plot
subplot(5,3,14);
plot(sc,Red_rec_N_correspondingtotalKL_norm,'-ys',sc,Red_rec_N2_correspondingtotalKL_norm,'-ms',sc,Red_rec_N3_correspondingtotalKL_norm,'-cs',sc,Red_rec_N4_correspondingtotalKL_norm,'-ks');
title('Red channel, corresponding condition of all non recorded neurons');
%Conversion rate 4 non recorded neurons plot
subplot(5,3,15);
plot(sc,Conversion_rec_N_correspondingtotalKL_norm,'-ys',sc,Conversion_rec_N2_correspondingtotalKL_norm,'-ms',sc,Conversion_rec_N3_correspondingtotalKL_norm,'-cs',sc,Conversion_rec_N4_correspondingtotalKL_norm,'-ks');
title('Conversion rate, corresponding condition of all non recorded neurons');



save (['normalized KL divergence' finalfile '.mat']);

end



