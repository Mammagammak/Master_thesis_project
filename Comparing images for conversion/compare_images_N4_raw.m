function compare_images_N4_raw(before,after,tags,finalfile)
%% What does the function do?
%This function enables user to compare fluroescence signal densitities between
%between images.
%Images are taken before any stimulation of the neurons in green and red channels
%these conditions are called before. Then the neuron is stimulated to fire predetermined 
%amount of action potentials and the after photo is taken for each condition.
%in which, before and after conditions for neurons allows the user to compare the signal intensity across
% 
%In this function also non-patched (the neurons that are not directly stimulated via electrode)
%are analyzed. This allows for the further analysis of connectivity within a neuronal circuit. 

%% INPUTS: 
% BEFORE and AFTER are structural arrays with at least one string, name of 
% the image files.  
% TAGS is a structural array with legends to be used for the images.
% FINALFILE is the name you would like to catalog the data under.  This is
% typically the name of the experiment. 

%4 neurons need to be chosen and entered for this function: 1 recorded/patched (stimulated with electrode)
% and 3 non-recorded/non-patched (not directly stimulated neurons)

%% Raw data used
%refers to the all pixel values across the experiment taken as their
%absolute values, without any processing
%
%% Sample entry
% compare_images_N3_raw ({'160215_E7_Preconversion_Green.tif';'160215_E7_Preconversion_Red.tif'},{'160215_E7_Postconversion_Green_Pulsed light.tif';'160215_E7_Postconversion_Red_Pulsed light.tif'},{'Green';'Red'},'160215')

%Cemre Kizilarmut (cemrekizilarmut@gmail.com)

close all


%% Choosing ROI, two choices will be made by the user: one for the recorded neuron one for the non recorded neuron. this ROI will be chosen to plot the data limited to taht ROI.
roi_recordedneuron_columns=input('enter the columns of the region of recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_recordedneuron_rows=input('enter the rows of the region of recorded neuron to analyze in an array. Exp:[y0 y1] = ');

roi_nonrecordedneuron_columns=input('enter the columns of the region of non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron_rows=input('enter the rows of the region of non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');

roi_nonrecordedneuron2_columns=input('enter the columns of the region of the second non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron2_rows=input('enter the rows of the region of the second non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');

roi_nonrecordedneuron3_columns=input('enter the columns of the region of the third non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron3_rows=input('enter the rows of the region of the third non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');


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
    for i=1:size(after,1);
    data.raw.afterrecorded{1,i}=data.raw.after{1,i}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    end
%first non recorded neuron roi extraction
    data.raw.beforenonrecorded{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    data.raw.beforenonrecorded{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
i=1;
    for i=1:size(after,1);
    data.raw.afternonrecorded{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    end
%second non recorded neuron roi extraction
    data.raw.beforenonrecorded2{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    data.raw.beforenonrecorded2{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.raw.afternonrecorded2{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    end
%third non recorded neuron roi extraction
    data.raw.beforenonrecorded3{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    data.raw.beforenonrecorded3{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.raw.afternonrecorded3{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    end

% recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
d=1;
b=1;
n=size(after,1)/2;
for k=1:n;
    G(1,b)=mean2(data.raw.afterrecorded{1,d});  
    d=d+2
    b=b+1
end




e=1;
f=2;
m=size(after,1)/2;    
for l=1:m;
    R(1,e)=mean2(data.raw.afterrecorded{1,f});
    e=e+1
    f=f+2
end


% non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi


g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
    G_N(1,g)=mean2(data.raw.afternonrecorded{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;    
for s=1:p;
    R_N(1,j)=mean2(data.raw.afternonrecorded{1,o});
    j=j+1
    o=o+2
end

% second non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
     G_N2(1,g)=mean2(data.raw.afternonrecorded2{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;   
for s=1:p;
    R_N2(1,j)=mean2(data.raw.afternonrecorded2{1,o});
    j=j+1
    o=o+2
end

% third non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
    G_N3(1,g)=mean2(data.raw.afternonrecorded3{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;    
for s=1:p;
    R_N3(1,j)=mean2(data.raw.afternonrecorded3{1,o});
    j=j+1
    o=o+2
end



if size(after,1)/2<9;
    v=1;
    y=9-(size(after,1)/2);
    for v=1:y;
    G(1,(size(after,1)/2+v))=0;
    R(1,(size(after,1)/2+v))=0;
    G_N(1,(size(after,1)/2+v))=0;
    R_N(1,(size(after,1)/2+v))=0;
    G_N2(1,(size(after,1)/2+v))=0;
    R_N2(1,(size(after,1)/2+v))=0;
    G_N3(1,(size(after,1)/2+v))=0;
    R_N3(1,(size(after,1)/2+v))=0;
    end
end
G_difference=(G)-(((G_N)+(G_N2)+(G_N3))/3);
R_difference=(R)-(((R_N)+(R_N2)+(R_N3))/3);

save (['Compare Images Time Series Raw Data for 4 neurons' finalfile '.mat']); 
x=[10 30 80 180 330 880 1530 2380 3430];
plot(x,G,'m-o',x,R,'b-o',x,G_N,'g--',x,R_N,'r--',x,G_difference,'y--',x,R_difference,'k--',x,G_N2,'g--',x,R_N2,'r--',x,G_N3,'g--',x,R_N3,'r--');
title('Raw Fluorescence level without any processing ');
xlabel('Cumulative Action Potentials Fired');
ylabel('Fluorescence Level');
legend('Green fluorescence recorded neuron','Red fluorescence recorded neuron','Green fluorescence non recorded neurons','Red fluorescence non recorded neurons','Green Fluorescence trend line','Red Fluorescence trend line');

