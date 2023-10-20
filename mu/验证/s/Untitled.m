gama1=2;
gama2=0.5;
for i=1:100
for j=1:100
h(i,j)=(1-(1-gama1)*PA1_ave(1,j)-(1-gama2)*PA2_ave(1,j));
end
  end