%% certain UV light power is needed in order to induce the conversion of
% of the green fluorescent protein. I have measured the light power corresponding 
% to each applied voltage over the LED. This allows me to predict and
% achieve the max light power without burning the LED
 
x=[18, 32, 52, 100, 148, 183, 200, 213];
y=[3, 6.3, 9.8, 17.4, 24.8, 26.6, 27, 29.4];
fivevolt=fittype('a+b*log(x)');
myfit=fit(x,y,fivevolt);
plot(myfit,x,y);
d=[18, 31, 32, 52, 66, 100, 122, 148, 183, 200, 213, 240, 270, 335, 360];
e=[3,7.1, 6.3, 9.8, 11.3, 17.4, 20.8, 24.8, 26.6, 27, 29.4, 27, 27, 31.5, 29];
adaptor=fittype('a+b*log(d)');
myfit2=fit(d',e',adaptor);
hold on;
plot(myfit2,d,e);
k=[31, 66, 122, 240, 270, 335, 360];
l=[7.1,11.3, 20.8, 27, 27, 31.5, 29];
together=fittype('a+b*log(k)');
myfit3=fit(k',l',together);
hold on;
plot(myfit3,k,l);

