function compare_images_N5 (before,after,tags,finalfile)
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
%
%Normalization is performed on the data by 
%% INPUTS: 
% BEFORE and AFTER are structural arrays with at least one string, name of 
% the image files.  
% TAGS is a structural array with legends to be used for the images.
% FINALFILE is the name you would like to catalog the data under.  This is
% typically the name of the experiment. 

%5 neurons need to be chosen and entered for this function: 1 recorded/patched (stimulated with electrode)
% and 4 non-recorded/non-patched (not directly stimulated neurons)

%% Data normalized by absolute max value
%Normalization was achieved by multiplying the values of each pixel by 2 and dividing each pixel by the 
%ultimate max value within that image (maximum within the same condition and the same channel)
%
%% Sample entry
% compare_images_N3 ({'160215_E7_Preconversion_Green.tif';'160215_E7_Preconversion_Red.tif'},{'160215_E7_Postconversion_Green_Pulsed light.tif';'160215_E7_Postconversion_Red_Pulsed light.tif'},{'Green';'Red'},'160215')

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

roi_nonrecordedneuron4_columns=input('enter the columns of the region of the fourth non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron4_rows=input('enter the rows of the region of the fourth non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
%% input control
    if size (before{1},1) ~= size (after{1},1)
     disp ('the number of images in the before and after conditions do not match!')
        return
    end

%% data import
    for lp=1:size(before,1)
        data.raw.before{lp} = imread(before{lp});
        data.raw.before{lp}=double(data.raw.before{lp});
        data.normal.before{lp} = data.raw.before{lp}./max(max(data.raw.before{lp}));
        
    end
    for a=1:size(after,1)
        data.raw.after{a} = imread(after{a});
        data.raw.after{a}=double(data.raw.after{a});
        data.normal.after{a} = data.raw.after{a}./max(max(data.raw.after{a}));
    end
    


data.labels = tags; 
data.inputs.before = before; 
data.inputs.after = after;



%% plot data in time 
%recorded neuron roi extraction   
    data.normal.beforerecorded{1,1}=data.normal.before{1,1}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    data.normal.beforerecorded{1,2}=data.normal.before{1,2}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    for i=1:size(after,1);
    data.normal.afterrecorded{1,i}=data.normal.after{1,i}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    end
%first non recorded neuron roi extraction
    data.normal.beforenonrecorded{1,1}=data.normal.before{1,1}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    data.normal.beforenonrecorded{1,2}=data.normal.before{1,2}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
i=1;
    for i=1:size(after,1);
    data.normal.afternonrecorded{1,i}=data.normal.after{1,i}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    end
%second non recorded neuron roi extraction
    data.normal.beforenonrecorded2{1,1}=data.normal.before{1,1}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    data.normal.beforenonrecorded2{1,2}=data.normal.before{1,2}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.normal.afternonrecorded2{1,i}=data.normal.after{1,i}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    end
%third non recorded neuron roi extraction
    data.normal.beforenonrecorded3{1,1}=data.normal.before{1,1}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    data.normal.beforenonrecorded3{1,2}=data.normal.before{1,2}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.normal.afternonrecorded3{1,i}=data.normal.after{1,i}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    end
%fourth non recorded neuron roi extraction
    
    data.normal.beforenonrecorded4{1,1}=data.normal.before{1,1}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
    data.normal.beforenonrecorded4{1,2}=data.normal.before{1,2}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.normal.afternonrecorded4{1,i}=data.normal.after{1,i}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
    end
    
% recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
d=1;
b=1;
n=size(after,1)/2;
for k=1:n;
    data.normal.afterrecorded{1,d}=(data.normal.afterrecorded{1,d}./data.normal.beforerecorded{1,1});
    G(1,b)=mean2(data.normal.afterrecorded{1,d});  
    d=d+2
    b=b+1
end




e=1;
f=2;
m=size(after,1)/2;    
for l=1:m;
    R(1,e)=mean2(data.normal.afterrecorded{1,f}./data.normal.beforerecorded{1,2});
    e=e+1
    f=f+2
end


% non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi


g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
    data.normal.afternonrecorded{1,h}=(data.normal.afternonrecorded{1,h}./data.normal.beforenonrecorded{1,1});
    G_N(1,g)=mean2(data.normal.afternonrecorded{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;    
for s=1:p;
    R_N(1,j)=mean2(data.normal.afternonrecorded{1,o}./data.normal.beforenonrecorded{1,2});
    j=j+1
    o=o+2
end

% second non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
    data.normal.afternonrecorded2{1,h}=(data.normal.afternonrecorded2{1,h}./data.normal.beforenonrecorded2{1,1});
    G_N2(1,g)=mean2(data.normal.afternonrecorded2{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;   
for s=1:p;
    R_N2(1,j)=mean2(data.normal.afternonrecorded2{1,o}./data.normal.beforenonrecorded2{1,2});
    j=j+1
    o=o+2
end

% third non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
    data.normal.afternonrecorded3{1,h}=(data.normal.afternonrecorded3{1,h}./data.normal.beforenonrecorded3{1,1});
    G_N3(1,g)=mean2(data.normal.afternonrecorded3{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;    
for s=1:p;
    R_N3(1,j)=mean2(data.normal.afternonrecorded3{1,o}./data.normal.beforenonrecorded3{1,2});
    j=j+1
    o=o+2
end

% fourth non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi

g=1;
h=1;
z=size(after,1)/2;
for t=1:z;
    data.normal.afternonrecorded4{1,h}=(data.normal.afternonrecorded4{1,h}./data.normal.beforenonrecorded4{1,1});
    G_N4(1,g)=mean2(data.normal.afternonrecorded4{1,h});  
    h=h+2
    g=g+1
end

j=1;
o=2;
p=size(after,1)/2;  
for s=1:p;
    R_N4(1,j)=mean2(data.normal.afternonrecorded4{1,o}./data.normal.beforenonrecorded4{1,2});
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
    G_N4(1,(size(after,1)/2+v))=0;
    R_N4(1,(size(after,1)/2+v))=0;
    end
end

G_difference=(G)-(((G_N)+(G_N2)+(G_N3)+(G_N4))/4);
R_difference=(R)-(((R_N)+(R_N2)+(R_N3)+(R_N4))/4);



save (['Compare Images Time Series Normalized for each condition for 5 neurons_' finalfile '.mat']); 
x=[10 30 80 180 330 880 1530 2380 3430];
plot(x,G,'m-o',x,R,'b-o',x,G_N,'g--',x,R_N,'r--',x,G_difference,'y--',x,R_difference,'k--',x,G_N2,'g--',x,R_N2,'r--',x,G_N3,'g--',x,R_N3,'r--',x,G_N4,'g--',x,R_N4,'r--');
title('After/Before conversion rate');
xlabel('Cumulative Action Potentials Fired');
ylabel('Conversion Rate');
legend('Green fluorescence recorded neuron','Red fluorescence recorded neuron','Green fluorescence non recorded neurons','Red fluorescence non recorded neurons','Green Fluorescence trend line','Red Fluorescence trend line');

