%This function is used to calculate how is the conversion from 
%green fluorescent protein to red fluorescent protein relates to 
%different parameters (size of the neuron, initial pixel value - which is the 
%count for fluorescent protein -, spikes induced to neurons, and duration of light exposure to the neuron)
function [lambda,psi,T,stats,F]=calculate_factor_analysis(factor_analysis);
directory=input('choose the directory=');

oldFolder = cd(directory);
s='.mat';
for_load=strcat(factor_analysis,s);
load(for_load);

%difference in the red channel of the neuron multiplied by the size of the
%neuron
m=size(Red_pre_diff_total_values_multiplied_size,2);
for i=1:m
difference{i}=Red_pre_diff_total_values_multiplied_size{i};
end

%assigning size of the neurons to a matrice 
size_neurons=b;
%assigning the normalized initial brightness of the neuron in the green channel
%to a matrice
initial_neurons=normalized_green_initial_values;

cd(oldFolder)
new_directory=input('enter the directory of without patch data you want to add=');

oldFolder=cd(new_directory);
to_load=strcat(new_directory,s);
load('conversion rates merged difference pixel wise no nucleus recorded neuron high sensitivity NB wo patch.mat');

n=size(Red_pre_diff_total_values_multiplied_size,2);
for i=1:n
difference_wo{i}=Red_pre_diff_total_values_multiplied_size{i};
end
size_neurons_wo=b;
initial_neurons_wo=normalized_green_initial_values;

% for i=1:size(b,2)
% K(i,:)=[b(i), normalized_green_initial_values(i), difference(i)];
% end
% m=1;
% 
% Y1=K(1:6,:);
% Y2=K(8:14,:);
% X=[Y1;Y2];
% 
% % bringing mean to zero
% mean_size=mean(X(:,1));
% mean_initial=mean(X(:,2));
% mean_difference=mean(X(:,3));
% Z(:,1)=X(:,1)-(mean_size);
% Z(:,2)=X(:,2)-(mean_initial);
% Z(:,3)=X(:,3)-(mean_difference);
% % calculating covariance matrices
% C1=cov(Z(:,1),Z(:,2));
% C2=cov(Z(:,2),Z(:,3));
% C3=cov(Z(:,1),Z(:,3));
%calculating eignvectors and eigenvalues of covariance matrices

%anova based on, size, initial value, spikes, and light duration

difference_all=[difference,difference_wo];
difference_vert=vertcat(difference_all{:});
difference_tr=difference_vert';
difference_vertical=difference_tr(:);
differece_horizontal=difference_vertical';

size_neurons_all=[size_neurons,size_neurons_wo];
size_neurons_all_tr=size_neurons_all';
size_vert=horzcat(size_neurons_all_tr,size_neurons_all_tr,size_neurons_all_tr,size_neurons_all_tr,size_neurons_all_tr,size_neurons_all_tr,size_neurons_all_tr,size_neurons_all_tr);
size_tr=size_vert';
size_neurons_horz=(size_tr(:))';

initial_neurons_all=[initial_neurons,initial_neurons_wo];
initial_neurons_all_tr=initial_neurons_all';
initial_vert=horzcat(initial_neurons_all_tr,initial_neurons_all_tr,initial_neurons_all_tr,initial_neurons_all_tr,initial_neurons_all_tr,initial_neurons_all_tr,initial_neurons_all_tr,initial_neurons_all_tr);
initial_tr=initial_vert';
initial_neurons_horz=(initial_tr(:)');

%spikes
s=[10 30 80 180 430 880 1530 2380];
spikes=[10 30 80 180 430 880 1530 2380];
for i=1:(m-1)
    spikes=[spikes,s];
end
ns=[0 0 0 0 0 0 0 0];
nspikes=[0 0 0 0 0 0 0 0];
for i=1:(n-1)
    nspikes=[nspikes,ns];
end
spikes_horz=[spikes,nspikes];

%light duration 
lights=[10 30 80 180 430 880 1530 2380];
nlights=[10 30 80 180 430 880 1530 2380];
for i=1:(m+n)-1
    nlights=[nlights,lights];
end



% y=[]
% 
% size=[]
% initial=[]
% light_duration=[]
% spike=[];
[p,tbl,stats,terms] = anovan(differece_horizontal,{size_neurons_horz,initial_neurons_horz,spikes_horz,nlights})

[lambda,psi,T,stats,F] = factoran(X,m);

[COEFF,SCORE,latentp,tsquare] = princomp(X)
[coeff,score,latent,tsquared,explained,mu] = pca(X)

save(['Factor analysis' factor_analysis '.mat']);
cd(oldFolder);
end
