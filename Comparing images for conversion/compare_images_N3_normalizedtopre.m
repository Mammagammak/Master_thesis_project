function compare_images_N3_normalizedtopre (before,after,tags,finalfile)
%% This function enables user to compare fluroescence signal densitities between
%between images.
%Images are taken before any stimulation of the neurons in green and red channels
%these conditions are called before. Then the neuron is stimulated to fire predetermined 
%amount of action potentials and the after photo is taken for each condition.
%in which, before and after conditions for neurons allows the user to compare the signal intensity across
% 
%In this function also non-patched (the neurons that are not directly stimulated via electrode)
%are analyzed. This allows for the further analysis of connectivity within a neuronal circuit. 
%
%% INPUTS: 
% BEFORE and AFTER are structural arrays with at least one string, name of 
% the image files.  
% TAGS is a structural array with legends to be used for the images.
% FINALFILE is the name you would like to catalog the data under.  This is
% typically the name of the experiment. 
%
%% Data normalized to previous condition is used
% refers to the all pixel values across the experiment are
%normalized to the before image's pixel values, in order to work with
%proportionated values.
%
%3 neurons need to be chosen and entered for this function: 1 recorded/patched (stimulated with electrode)
% and 2 non-recorded/non-patched (not directly stimulated neurons)

%% Sample entry
% compare_images_N3_normalizedtopre ({'160215_E7_Preconversion_Green.tif';'160215_E7_Preconversion_Red.tif'},{'160215_E7_Postconversion_Green_Pulsed light.tif';'160215_E7_Postconversion_Red_Pulsed light.tif'},{'Green';'Red'},'160215')

%Cemre Kizilarmut (cemrekizilarmut@gmail.com)

close all

%% Choosing ROI, two choices will be made by the user: one for the recorded neuron one for the non recorded neuron. this ROI will be chosen to plot the data limited to taht ROI.
roi_recordedneuron_columns=input('enter the columns of the region of recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_recordedneuron_rows=input('enter the rows of the region of recorded neuron to analyze in an array. Exp:[y0 y1] = ');

roi_nonrecordedneuron_columns=input('enter the columns of the region of non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron_rows=input('enter the rows of the region of non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');

roi_nonrecordedneuron2_columns=input('enter the columns of the region of the second non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron2_rows=input('enter the rows of the region of the second non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');


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
        b=mod(a,2);
        if  b==0;
            b
        data.normal.after{a} = data.raw.after{a}./max(max(data.raw.before{2}));
        else
           data.normal.after{a} = data.raw.after{a}./max(max(data.raw.before{1})); 
        end
        
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


    
% first condition to be plotted in the graph, pre condition without any electrical stimulation    
pregreen=mean2(data.normal.beforerecorded{1,1}./data.normal.beforerecorded{1,1});
prered=mean2(data.normal.beforerecorded{1,2}./data.normal.beforerecorded{1,2});    
    
% recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi
G(1,1)=pregreen;
d=1;
b=2;
n=size(after,1)/2;
for k=1:n;
    data.normal.afterrecorded{1,d}=(data.normal.afterrecorded{1,d}./data.normal.beforerecorded{1,1});
    G(1,b)=mean2(data.normal.afterrecorded{1,d});  
    d=d+2
    b=b+1
end



R(1,1)=prered;
e=2;
f=2;
m=size(after,1)/2;    
for l=1:m;
    R(1,e)=mean2(data.normal.afterrecorded{1,f}./data.normal.beforerecorded{1,2});
    e=e+1
    f=f+2
end


% non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi

G_N(1,1)=pregreen;
g=2;
h=1;
z=size(after,1)/2;
for t=1:z;
    data.normal.afternonrecorded{1,h}=(data.normal.afternonrecorded{1,h}./data.normal.beforenonrecorded{1,1});
    G_N(1,g)=mean2(data.normal.afternonrecorded{1,h});  
    h=h+2
    g=g+1
end

R_N(1,1)=prered;

j=2;
o=2;
p=size(after,1)/2;    
for s=1:p;
    R_N(1,j)=mean2(data.normal.afternonrecorded{1,o}./data.normal.beforenonrecorded{1,2});
    j=j+1
    o=o+2
end

% second non recorded neuron after/before calculation separately for green and red, and later mean is taken as representation of the roi

G_N2(1,1)=pregreen;
g=2;
h=1;
z=size(after,1)/2;
for t=1:z;
    data.normal.afternonrecorded2{1,h}=(data.normal.afternonrecorded2{1,h}./data.normal.beforenonrecorded2{1,1});
    G_N2(1,g)=mean2(data.normal.afternonrecorded2{1,h});  
    h=h+2
    g=g+1
end

R_N2(1,1)=prered;
j=2;
o=2;
p=size(after,1)/2;   
for s=1:p;
    R_N2(1,j)=mean2(data.normal.afternonrecorded2{1,o}./data.normal.beforenonrecorded2{1,2});
    j=j+1
    o=o+2
end

% To observe the ratiometric property of the CaMPARI, totalfluorecence in
% each condition is calculated by summing up the red and green fluorescence level at that condition. 

xpre=[0 10 30 80 180 330 880 1530 2380 3430];
totalfluorescence(1,1)=mean2(data.normal.beforerecorded{1,1})+mean2(data.normal.beforerecorded{1,2})
wz=2;
for w=1:2:size(after,1);
totalfluorescence(1,wz)=mean2(data.normal.afterrecorded{1,w})+mean2(data.normal.afterrecorded{1,(w+1)})
wz=wz+1;
end
totalfluorescence

% In the case of an experiment being incomplete this part adds zero's to
% the end of the matrices to make sure the lengths are same.



if size(after,1)/2<9;
    v=2;
    y=10-(size(after,1)/2);
    for v=1:y;
    G(1,(size(after,1)/2+v))=0;
    R(1,(size(after,1)/2+v))=0;
    G_N(1,(size(after,1)/2+v))=0;
    R_N(1,(size(after,1)/2+v))=0;
    G_N2(1,(size(after,1)/2+v))=0;
    R_N2(1,(size(after,1)/2+v))=0;
    totalfluorescence(1,(size(after,1)/2+v))=0;
    end
end


% difference matrices give the net fluorescence change between the recorded
% neuron and the non-recorded neurons
G_difference=(G)-(((G_N)+(G_N2))/2);
R_difference=(R)-(((R_N)+(R_N2))/2);

% Gradient is calculated in order to determine how different the recorded
% neurons is from the non-recorded neurons. The values are saved to a csv
% file for later to be analyzed statistically

%csv file need to be created by the matlab prior to starting registering
%real values. if the file is created manually with dummy values matlab
%gives error while reading it. so the first row is always to be disgarded.
[Frecorded]=gradient(R);
[Fnonrecorded]=gradient(((R_N)+(R_N2))/2);
slope=max(Frecorded);

filename1=('Frecorded and Fnonrecorded.csv');
M=csvread(filename1);
M
Mfinal=[M;Frecorded;Fnonrecorded];
csvwrite(filename1,Mfinal);
Mfinal

% gradient values of the net red fluorescence (recorded neuron fluorescence
% level- non-recorded neurons averaged fluorescence level) exported into
% csv file.
[Frecordednet]=gradient(R_difference);

filename2=('Frecorded net.csv');
S=csvread(filename2);
S
Sfinal=[S;Frecordednet];
csvwrite(filename2,Sfinal);
Sfinal

% after/before conversion rate of recorded neuron and non recorded
% neurons(averaged), registered to csv file.

[Frecordedvalue]=R;
[Fnonrecordedvalue]=(((R_N)+(R_N2))/2);

filename3=('Frecorded and Fnonrecorded absolute values.csv');
N=csvread(filename3);
N
Nfinal=[N;Frecordedvalue;Fnonrecordedvalue];
csvwrite(filename3,Nfinal);
Nfinal


save (['Compare Images time series normalized to pre condition for 3 neurons_' finalfile '.mat']); 
x=[0 10 30 80 180 330 880 1530 2380 3430];
plot(x,G,'m-o',x,R,'b-o',x,G_N,'g--',x,R_N,'r--',x,G_difference,'y--',x,R_difference,'k--',xpre,totalfluorescence,'c-*',x,G_N2,'g--',x,R_N2,'r--');
title('After/Before conversion rate' );
xlabel('Cumulative Action Potentials Fired');
ylabel('Conversion Rate');
legend('Green fluorescence recorded neuron','Red fluorescence recorded neuron','Green fluorescence non recorded neurons','Red fluorescence non recorded neurons','Green Fluorescence trend line','Red Fluorescence trend line');

