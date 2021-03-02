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
    Gmax=max(G);
b=1;
n=(size(after,1)/2)+1;
    for k=1:n;
    Gprime(1,b)=G(1,b)./Gmax;
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
    Rmax=max(R);
e=1;
m=(size(after,1)/2)+1;
    for l=1:m;
        Rprime(1,e)=R(1,e)./Rmax;
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
    G_Nmax=max(G_N);
b=1;
n=(size(after,1)/2)+1;
    for k=1:n;
    G_Nprime(1,b)=G_N(1,b)./G_Nmax;
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
    R_Nmax=max(R_N);
e=1;
m=(size(after,1)/2)+1;
    for l=1:m;
        R_Nprime(1,e)=R_N(1,e)./R_Nmax;
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
    G_N2max=max(G_N2);
b=1;
n=(size(after,1)/2)+1;
    for k=1:n;
    G_N2prime(1,b)=G_N2(1,b)./G_N2max;
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
    R_N2max=max(R_N2);
e=1;
m=(size(after,1)/2)+1;
    for l=1:m;
        R_N2prime(1,e)=R_N(1,e)./R_N2max;
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
    G_N3max=max(G_N3);
b=1;
n=(size(after,1)/2)+1;
    for k=1:n;
    G_N3prime(1,b)=G_N3(1,b)./G_N3max;
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
    R_N3max=max(R_N3);
e=1;
m=(size(after,1)/2)+1;
    for l=1:m;
        R_N3prime(1,e)=R_N3(1,e)./R_N3max;
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
    G_N4max=max(G_N4);
b=1;
n=(size(after,1)/2)+1;
    for k=1:n;
    G_N4prime(1,b)=G_N4(1,b)./G_N4max;
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
    R_N4max=max(R_N4);
e=1;
m=(size(after,1)/2)+1;
    for l=1:m;
        R_N4prime(1,e)=R_N4(1,e)./R_N4max;
        e=e+1;
    end


%% Rough formula for the behavior of Campari
%independent from [Ca+2] concentration. k is the contribution of the spike frequency
%to conversion, A is the effect each spike has to the fluorescence level.

%there are 10 conditions at the protocol. First condition is before any
%stimulation and it should give the baseline fluorescence level (Also the fluorescence level of that cell in baseline calcium concentration).

F=a.C(t)+beta+error
C(t)=C(t-1).*(1-(freqsamp/decaytimeCa))+(freqsamp/decaytimeCa).*Cb+A.*s





freq=input('enter the frequency of the APs fired(Hz)=');
s=[0 10 20 50 100 250 450 650 850 1050];
syms A;
syms k;

F=R;
Fb=F(1)
Fr(1)=0;


for i=2:10;
    Fr(i)=A.*s(i)+freq.*k+Fb-F(i);
end

Fr(1)==0;
Fr(2)==0;
Fr(3)==0;
Fr(4)==0;
Fr(5)==0;
Fr(6)==0;
Fr(7)==0;
Fr(8)==0;
Fr(9)==0;
Fr(10)==0;
    
    
%% Calculation of the fractional change in the signal intensity
Cpre=Rprime./Gprime;
Cpost=R_Nprime./G_Nprime;
Cpost2=R_N2prime./G_N2prime;
Cpost3=R_N3prime./G_N3prime;
Cpost4=R_N4prime./G_N4prime;

x=[0 10 30 80 180 330 880 1530 2380 3430];


%% Alternative calculation: Subtractional normalization's ratio to total fluorescence 

Cpres=(Rprime-Gprime)./(Rprime+Gprime);
Cposts=(R_Nprime-G_Nprime)./(R_Nprime+G_Nprime);
Cposts2=(R_N2prime-G_N2prime)./(R_N2prime+G_N2prime);
Cposts3=(R_N3prime-G_N3prime)./(R_N3prime+G_N3prime);
Cposts4=(R_N4prime-G_N4prime)./(R_N4prime+G_N4prime);

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


