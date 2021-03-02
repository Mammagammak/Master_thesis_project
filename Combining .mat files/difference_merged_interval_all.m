function difference_merged_interval_all(diff_values)
% this function merges the interested *.mat files for further analysis
directory_to_save=input('enter the name of the folder you wish to save the current file to: e.g. "Conversion rates medium sensitivity NB w patch merged". = ');
oldFolder = cd(diff_values);
interval_spike_counts=input('Enter the spike counts of the interval varying experiment you want to analyze as a numeric variable: (if the experiment has standard protcol enter 0)=');

cd(interval_spike_counts);

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

all_ROI{i}=m.ROI; % ROI from each experiment, this matrice contains 1's as the ROI and 0's as the other regions in the image
a{i}=find(all_ROI{i}==1); %extract the ROI 
b(i)=size(a{i},1); %size of the neuron as pixel count from each experiment

all_diff_total_values{i}=m.all_diff_values; %all conditions of recorded neuron from each experiment

Red_pre_diff_total_values{i}=all_diff_total_values{i}{3}; % red pre condition conversion values from each experiment

clearvars m;
end    
    

 
for j=1:size(Red_pre_diff_total_values,2)
    carrier=Red_pre_diff_total_values{j};
    carrier=[0,carrier];
    Red_pre_diff_total_values_zero_added{j}=carrier;
end   

s100=[0 100 200 300 400 500];
s120=[0 120 240 360 480 600];
s135=[0 135 270 405 540 675];
s150=[0 150 300 450 600 750];
s175=[0 175 350 525 700 875];
s250=[0 250 500 750 1000 1250];
s450=[0 450 900 1350 1800 2250];
X=[s100,s100,s250,s250,s250,s250,s250,s450,s450,s175,s150,s120,s135,s175,s450,s175,s175,s250,s450,s450,s175,s175,s450,s450,s450];



spike=0;
Y=[];

f1=figure;
for i=1:size(Red_pre_diff_total_values_zero_added,2)
     y_carrier=Red_pre_diff_total_values_zero_added{i};
     Y=[Y,y_carrier];
end


figure(f1) 
h=stem(X,Y,'LineStyle','none','MarkerFaceColor','black', 'MarkerEdgeColor','black');
hold on 



for i=1:size(Red_pre_diff_total_values_zero_added,2)
     y_carrier=Red_pre_diff_total_values_zero_added{i};
     
if i==1 %100 red
    %h(i).Color='red' ;
    Color='Red' ;
    spike=s100;
elseif i==2 
        %h(i).Color='red';
        Color='Red' ;
        spike=s100;
        
elseif i==3 % 250
        %h(i).Color= 'black';
        Color= 'k';
        spike=s250;
elseif i==4
        %h(i).Color= 'black';
        Color= 'k';
        spike=s250;
elseif i==5
        %h(i).Color='black' ;
        Color= 'k';
        spike=s250;
elseif i==6
        %h(i).Color='black';
        Color= 'k';
        spike=s250;
elseif i==7
         %h(i).Color='black';
         Color= 'k';
          spike=s250;
elseif i==8 % 450
    %h(i).Color='magenta';
    Color='m';
    spike=s450;
    elseif i==9
            %h(i).Color='magenta';
            Color='m';
            spike=s450;
            elseif i==10 % 175 blue
              % h(i).Color='blue';
               Color='b';
                spike=s175;
                elseif i==11 % 150 cyan
                   % h(i).Color='cyan';
                    Color='c';
                    spike=s150;
                    elseif i==12 %120 yellow
                    %    h(i).Color='yellow';
                        Color='y';
                        spike=s120;
                        elseif i==13 % 135 green
                       %     h(i).Color='green';
                            Color='g';
                            spike=s135;
                            elseif i==14
                          %      h(i).Color='blue';
                                Color='b';
                                spike=s175;
                                elseif i==15
                            %    h(i).Color='magenta';
                                Color='m';
                                spike=s450;
                                elseif i==16
                             %   h(i).Color='blue';
                                Color='b';
                                spike=s175;
                                elseif i==17
                            %    h(i).Color='blue';
                                Color='b';
                                spike=s175;
                                elseif i==18
                             %   h(i).Color='black';
                                Color='k';
                                spike=s250;
                                elseif i==19
                             %  h(i).Color='magenta';
                               Color='m';
                                spike=s450;
                                 elseif i==20
                               % h(i).Color='magenta';
                                Color='m';
                                spike=s450;
                                 elseif i==21
                             %   h(i).Color='blue';
                                Color='b';
                                spike=s175;
                                 elseif i==22
                              % h(i).Color='blue';
                               Color='b';
                                spike=s175;
                                 elseif i==23
                              %  h(i).Color='magenta';
                                Color='m';
                                spike=s450;
                                 elseif i==24
                              % h(i).Color='magenta';
                               Color='m';
                                spike=s450;
                                 elseif i==25
                               % h(i).Color='magenta';
                                Color='m';
                                spike=s450;
end

plot(spike,y_carrier,Color);
hold on



end

hold off



%% saving variables



filename1='interval varying all spikes';


print(f1,filename1,'-dsvg','-r600');


save([filename1 '.mat']);
cd(oldFolder);
end

