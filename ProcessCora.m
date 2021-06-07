clear;
clc;

load('Cora.mat');


n = size(Label,1);
X= zeros(n,1);
for i =1:n

X(i) = find(Label(i,:)>0);
end
Label = X;
Attributes = double(Attributes);

save('Cora.mat','Network', 'Attributes', 'Label');