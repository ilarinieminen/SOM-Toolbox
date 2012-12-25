

clf reset;
echo on

% This script demonstrates the use of the RSOM.
clc;

% load the example dissimilarity data
load exampleDissimilarity.mat;

% display the eigenvalue spectrum
[V, D] = eig(Dissim);
eVals = diag(D);
figure; bar(eVals);


pause % Strike any key to continues...
clc;


% init RSOM
sMap = rsom_lininit(Dissim,  [10 10]);

% train the RSOM
sMap = rsom_batchtrain(sMap, Dissim);


% Display the U-Matrix
figure;
rsom_show(sMap, Dissim);

% The U-Matrix shows clearly the structure of the data, i.e. that two
% clusters are available

pause % Strike any key to continues...
clc;

% do linear embedding of the distance matrix into a 3 dimensional space
x = cmdscale(Dissim.^(1/2));
x = x(:,1:3);

% plot the resulting data
h = figure; hold on;
plot3(x(1:100, 1), x(1:100, 2), x(1:100, 3), '.b');
plot3(x(101:200, 1), x(101:200, 2), x(101:200, 3), '.r');

% Since approximated vectorial data are available now, we can compute
% the (approximated) neuron positions and plot them
Neurons = sMap.cCodebook * x;
figure(h);

% create a som struct
sMapSOM = som_map_struct(3, sMap.topol); 
som_grid(sMapSOM,'Coord',Neurons);

echo off;


