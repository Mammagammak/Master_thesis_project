% Smoothing spline fit type is used for the analysis of the KL divergence
% data 
z=[10 20 50 100 250 450 650 850]
c=[0.0105526263139047,0.000951939841187198,0.000157905413659733,0.000293702307496319,0.127118382687927,8.18080207453548,1.48133635285251,2.20024502447044];
a=[0.0105526263139047,0.0162132003675911,0.0164702649844923,0.0175525350377225,0.232237392240152,10.2138600049918,17.6324307211556,21.0309991831404];
[xData, yData] = prepareCurveData( z, a );
t=(10:1:850);
coloumnvector=vertcat(t);
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.95;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


% syms f(x)
% f(x)= x.^3.*coeffs(1,1)+x.^2.*coeffs(1,2)+x.*coeffs(1,3)+coeffs(1,4);

[fx, fxx]=differentiate(fitresult,coloumnvector);
feval(fitresult,xData)
fastest_change_y_value=max(fx);
idx_max=(fx==fastest_change_y_value)
fastest_change_spike_value=coloumnvector(idx_max);
% fit_formula=poly2sym(w,k);
% ezplot(f(x), [10 19]);
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
subplot(3,1,1)
plot( fitresult, xData, yData );
subplot(3,1,2)
plot(t,fx)
subplot(3,1,3)
plot(t,fxx)



%%determining the dynamic range 
% the largest range in which the first derivative stays postive (the widest range that the original function is increasing)
% there are experiments where the conversion starts early and even though
% as a result of fit the dfunction seems decreasing during intervals, but
% this could be disgarded because small spike conversions are not
% consistent througout all experiments.
sign_fx=sign(fx);
idx_positive=(sign_fx==1);
%extracting the biggest uninterrupted positive region
o=1;
s=1;
for i=1:size(idx_positive,1)
    if i<size(idx_positive,1) & idx_positive((i+1),1)~=idx_positive(i,1) & idx_positive((i+1),1)==1
            change_from_zero_to_one{o}=i;
            o=o+1;
    elseif i<size(idx_positive,1) & idx_positive((i+1),1)~=idx_positive(i,1) & idx_positive((i+1),1)==0
            change_from_one_to_zero{s}=i;
            s=s+1;
    end
end

last_zero_to_one_change=max(change_from_zero_to_one{:});
last_one_to_zero_change=max(change_from_one_to_zero{:});
if last_zero_to_one_change>last_one_to_zero_change
    dynamic_range=[last_zero_to_one_change,850];
else
    dynamic_range=[last_zero_to_one_change,last_one_to_zero_change];
end


