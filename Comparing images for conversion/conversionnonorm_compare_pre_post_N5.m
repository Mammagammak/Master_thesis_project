function conversion_compare_pre_post_N5(before,after,tags,finalfile)
% compares the conversion rate of the presynaptic neuron with the
% postsynaptic neurons.

close all

%% Choosing ROI, coordinates will be entered by the user individually for each neuron; one for the recorded neuron and four for the non recorded neurons. This ROI will be used to limit the data process to that specific ROI.
roi_recordedneuron_columns=input('enter the columns of the region of recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_recordedneuron_rows=input('enter the rows of the region of recorded neuron to analyze in an array. Exp:[y0 y1] = ');
recorded_distance=0;

roi_nonrecordedneuron_columns=input('enter the columns of the region of non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron_rows=input('enter the rows of the region of non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded_distance=('enter the distance of the first non recorded neuron from the recorded neuron=');


roi_nonrecordedneuron2_columns=input('enter the columns of the region of the second non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron2_rows=input('enter the rows of the region of the second non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded2_distance=('enter the distance of the second non recorded neuron from the recorded neuron=');


roi_nonrecordedneuron3_columns=input('enter the columns of the region of the third non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron3_rows=input('enter the rows of the region of the third non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded3_distance=('enter the distance of the third non recorded neuron from the recorded neuron=');


roi_nonrecordedneuron4_columns=input('enter the columns of the region of the fourth non-recorded neuron to analyze in an array. Exp:[x0 x1] = ');
roi_nonrecordedneuron4_rows=input('enter the rows of the region of the fourth non-recorded neuron to analyze in an array. Exp:[y0 y1] = ');
nonrecorded4_distance=('enter the distance of the fourth non recorded neuron from the recorded neuron=');



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
%recorded neuron (presynaptic neuron) roi extraction (from this point all the process is made on these ROI's)   
    data.raw.beforerecorded{1,1}=data.raw.before{1,1}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    data.raw.beforerecorded{1,2}=data.raw.before{1,2}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    for i=1:size(after,1);
    data.raw.afterrecorded{1,i}=data.raw.after{1,i}(roi_recordedneuron_rows(1,1):roi_recordedneuron_rows(1,2),roi_recordedneuron_columns(1,1):roi_recordedneuron_columns(1,2));
    end
%first non recorded neuron (first post synaptic neuron) roi extraction
    data.raw.beforenonrecorded{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    data.raw.beforenonrecorded{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
i=1;
    for i=1:size(after,1);
    data.raw.afternonrecorded{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron_rows(1,1):roi_nonrecordedneuron_rows(1,2),roi_nonrecordedneuron_columns(1,1):roi_nonrecordedneuron_columns(1,2));
    end
%second non recorded neuron (second post synaptic neuron) roi extraction
    data.raw.beforenonrecorded2{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    data.raw.beforenonrecorded2{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.raw.afternonrecorded2{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron2_rows(1,1):roi_nonrecordedneuron2_rows(1,2),roi_nonrecordedneuron2_columns(1,1):roi_nonrecordedneuron2_columns(1,2));
    end
%third non recorded neuron (third post synaptic neuron) roi extraction
    data.raw.beforenonrecorded3{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    data.raw.beforenonrecorded3{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.raw.afternonrecorded3{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron3_rows(1,1):roi_nonrecordedneuron3_rows(1,2),roi_nonrecordedneuron3_columns(1,1):roi_nonrecordedneuron3_columns(1,2));
    end
%fourth non recorded neuron (fourth post synaptic neuron)roi extraction
    
    data.raw.beforenonrecorded4{1,1}=data.raw.before{1,1}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
    data.raw.beforenonrecorded4{1,2}=data.raw.before{1,2}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
i=1;   
    for i=1:size(after,1);
    data.raw.afternonrecorded4{1,i}=data.raw.after{1,i}(roi_nonrecordedneuron4_rows(1,1):roi_nonrecordedneuron4_rows(1,2),roi_nonrecordedneuron4_columns(1,1):roi_nonrecordedneuron4_columns(1,2));
    end
    
    
%% For each ROI, max pixel value among the whole time series is found. The whole values are divided by this max value for the channel normalization.

% Recorded neuron 
   
%pre condition assigned to the first value at the matrix   
G(1,1)=mean2(data.raw.beforerecorded{1,1});

d=1;
b=2;
n=size(after,1)/2;
    for k=1:n;
    G(1,b)=mean2(data.raw.afterrecorded{1,d}); 
     d=d+2;
     b=b+1;
    end
   

    
R(1,1)=mean2(data.raw.beforerecorded{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m;
    R(1,e)=mean2(data.raw.afterrecorded{1,f});
    f=f+2;
    e=e+1;
    end
   
           

% first non recorded neuron 

G_N(1,1)=mean2(data.raw.beforenonrecorded{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n;
    G_N(1,b)=mean2(data.raw.afternonrecorded{1,d}); 
     d=d+2;
     b=b+1;
    end
   
    
R_N(1,1)=mean2(data.raw.beforenonrecorded{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m;
    R_N(1,e)=mean2(data.raw.afternonrecorded{1,f});
    f=f+2;
    e=e+1;
    end
   


% second non recorded neuron 

G_N2(1,1)=mean2(data.raw.beforenonrecorded2{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n;
    G_N2(1,b)=mean2(data.raw.afternonrecorded2{1,d}); 
     d=d+2;
     b=b+1;
    end
    
    
R_N2(1,1)=mean2(data.raw.beforenonrecorded2{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m;
    R_N2(1,e)=mean2(data.raw.afternonrecorded2{1,f});
    f=f+2;
    e=e+1;
    end
   
% third non recorded neuron 

G_N3(1,1)=mean2(data.raw.beforenonrecorded3{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n;
    G_N3(1,b)=mean2(data.raw.afternonrecorded3{1,d}); 
     d=d+2;
     b=b+1;
    end
  
    
R_N3(1,1)=mean2(data.raw.beforenonrecorded3{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m;
    R_N3(1,e)=mean2(data.raw.afternonrecorded3{1,f});
    f=f+2;
    e=e+1;
    end
   
% fourth non recorded neuron 

G_N4(1,1)=mean2(data.raw.beforenonrecorded4{1,1});
d=1;
b=2;
n=size(after,1)/2;
    for k=1:n;
    G_N4(1,b)=mean2(data.raw.afternonrecorded4{1,d}); 
     d=d+2;
     b=b+1;
    end
    
    
R_N4(1,1)=mean2(data.raw.beforenonrecorded4{1,2});
e=2;
f=2;
m=size(after,1)/2;    
    for l=1:m;
    R_N4(1,e)=mean2(data.raw.afternonrecorded4{1,f});
    f=f+2;
    e=e+1;
    end
   

    
%% Calculation of the fractional change in the signal intensity
Cpre=R./G;
Cpost=R_N./G_N;
Cpost2=R_N2./G_N2;
Cpost3=R_N3./G_N3;
Cpost4=R_N4./G_N4;

x=[0 10 30 80 180 330 880 1530 2380 3430];


%% Alternative calculation: Subtractional normalization's ratio to total fluorescence 

Cpres=(R-G)./(R+G);
Cposts=(R_N-G_N)./(R_N+G_N);
Cposts2=(R_N2-G_N2)./(R_N2+G_N2);
Cposts3=(R_N3-G_N3)./(R_N3+G_N3);
Cposts4=(R_N4-G_N4)./(R_N4+G_N4);

%% Print the value of the distances

recorded_distance
nonrecorded_distance
nonrecorded2_distance
nonrecorded3_distance
nonrecorded4_distance



%% Plot and Save

save (['Conversion compare between presynaptic and postsynaptic neurons no background normalization_5 neurons_' finalfile '.mat']); 
x=[0 10 30 80 180 330 880 1530 2380 3430];
plot(x,Cpre,'m-o',x,Cpost,'b-o',x,Cpost2,'g-o',x,Cpost3,'r-o',x,Cpost4,'y-o',x,Cpres,'m-*',x,Cposts,'b-*',x,Cposts2,'g-*',x,Cposts3,'r-*',x,Cposts4,'y-*');

title('Conversion rate vs. spikes for differently spaced neurons' );
xlabel('Cumulative Action Potentials Fired');
ylabel('Conversion Rate');
legend('Recorded (pre) neuron ratiometric change of conversion (fractional)','NonRecorded (post)first neuron ratiometric change of conversion (fractional)','NonRecorded (post) second neuron ratiometric change of conversion (fractional)','NonRecorded (post) third neuron ratiometric change of conversion (fractional)','NonRecorded (post) fourth neuron ratiometric change of conversion (fractional)','Recorded (pre) neuron ratiometric change of conversion (subtractional)','NonRecorded (post) first neuron ratiometric change of conversion (subtractional)','NonRecorded (post) second neuron ratiometric change of conversion (subtractional)','NonRecorded (post) third neuron ratiometric change of conversion (subtractional)','NonRecorded (post) fourth neuron ratiometric change of conversion (subtractional)');


