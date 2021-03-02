%testing distance contribution on KL divergence values
%closest neurons are chosen as the neurons in the range of (8-50.48)
B=P<=50.48;
B

s = 1;
nB    = numel(B);
Match = cell(1, nB)
i=1;
for iB = 1:nB
   
  Match{iB} = find(B(iB) == s);
  if B(iB)==s;
     A{i}=iB;
     i=i+1;
  end
end
A
g=size(A,2);
l=1;
for l=1:g;
   c= mod(A(l),4);
   if c==0;
       exp_number=A(l)/4;
       
       if c==1;
           exp_number=(A(l)+3)/4;
           if c==2;
               exp_number=(A(l)+2)/4;
               if c==3;
                    exp_number=(A(l)+3)/4;
               end
           end
       end
   end
   C{l}=exp_number;
   D{l}=c;
end

%extracting the interested rows based on the exp_number info
l=1;
j=size(C,2)
for l=1:j;
    h=C{l};
    if D{l}==0;
        W=((h-1).*24)+12;
        closest_neurons_red_pre_difference_KL{l}=T(W,:);
    else
        W=((h-1).*24)+8+D{l};
        closest_neurons_red_pre_difference_KL{l}=T(W,:);
    end
end

%furthest neurons
V=P>=492.72;
V

s = 1;
nB    = numel(V);
Match = cell(1, nV)
i=1;
for iV = 1:nV
   
  Match{iV} = find(V(iV) == s);
  if V(iV)==s;
     Q{i}=iV;
     i=i+1;
  end
end
Q
g=size(Q,2);
l=1;
for l=1:g;
   c= mod(Q(l),4);
   if c==0;
       exp_number=Q(l)/4;
       
       if c==1;
           exp_number=(Q(l)+3)/4;
           if c==2;
               exp_number=(Q(l)+2)/4;
               if c==3;
                    exp_number=(Q(l)+3)/4;
               end
           end
       end
   end
   C_F{l}=exp_number;
   D_F{l}=c;
end

%extracting the interested rows based on the exp_number info
l=1;
j=size(C,2)
for l=1:j;
    h=C_F{l};
    if D_F{l}==0;
       Y=((h-1).*24)+12;
        closest_neurons_red_pre_difference_KL{l}=T(Y,:);
    else
        Y=((h-1).*24)+8+D_F{l};
        furthest_neurons_red_pre_difference_KL{l}=T(Y,:);
    end
end
x=[10 20 50 100 250 450 650 850 1050];    
plot(x,closest_neurons_red_pre_difference_KL,'r',x,furthest_neurons_red_pre_difference_KL,'g')

