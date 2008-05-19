
%SOM_DEMO2 Basic usage of the SOM Toolbox.

% Contributed to SOM Toolbox 2.0, February 11th, 2000 by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/

% Version 1.0beta juuso 071197 
% Version 2.0beta juuso 070200

clf reset;
figure(gcf)
echo on



clc
%    ==========================================================
%    SOM_DEMO2 - BASIC USAGE OF SOM TOOLBOX
%    ==========================================================

%    som_data_struct    - Create a data struct.
%    som_read_data      - Read data from file.
%
%    som_normalize      - Normalize data.
%    som_denormalize    - Denormalize data.
%
%    som_make           - Initialize and train the map. 
%
%    som_show           - Visualize map.
%    som_show_add       - Add markers on som_show visualization.
%    som_grid           - Visualization with free coordinates.
%
%    som_autolabel      - Give labels to map.
%    som_hits           - Calculate hit histogram for the map.

%    BASIC USAGE OF THE SOM TOOLBOX

%    The basic usage of the SOM Toolbox proceeds like this: 
%      1. construct data set
%      2. normalize it
%      3. train the map
%      4. visualize map
%      5. analyse results

%    The four first items are - if default options are used - very
%    simple operations, each executable with a single command.  For
%    the last, several different kinds of functions are provided in
%    the Toolbox, but as the needs of analysis vary, a general default
%    function or procedure does not exist. 

pause % Strike any key to construct data...



clc
%    STEP 1: CONSTRUCT DATA
%    ======================

%    The SOM Toolbox has a special struct, called data struct, which
%    is used to group information regarding the data set in one
%    place.

%    Here, a data struct is created using function SOM_DATA_STRUCT.
%    First argument is the data matrix itself, then is the name 
%    given to the data set, and the names of the components
%    (variables) in the data matrix.

D = rand(1000,3); % 1000 samples from unit cube
sData = som_data_struct(D,'name','unit cube','comp_names',{'x','y','z'});
			
%    Another option is to read the data directly from an ASCII file.
%    Here, the IRIS data set is loaded from a file (please make sure
%    the file can be found from the current path):

try, 
  sDiris = som_read_data('iris.data');
catch
  echo off

  warning('File ''iris.data'' not found. Using simulated data instead.')
  
  D = randn(50,4); 
  D(:,1) = D(:,1)+5;     D(:,2) = D(:,2)+3.5; 
  D(:,3) = D(:,3)/2+1.5; D(:,4) = D(:,4)/2+0.3;
  D(find(D(:)<=0)) = 0.01; 
  
  D2 = randn(100,4); D2(:,2) = sort(D2(:,2));
  D2(:,1) = D2(:,1)+6.5; D2(:,2) = D2(:,2)+2.8; 
  D2(:,3) = D2(:,3)+5;   D2(:,4) = D2(:,4)/2+1.5;
  D2(find(D2(:)<=0)) = 0.01; 
  
  sDiris = som_data_struct([D; D2],'name','iris (simulated)',...
			  'comp_names',{'SepalL','SepalW','PetalL','PetalW'});
  sDiris = som_label(sDiris,'add',[1:50]','Setosa');
  sDiris = som_label(sDiris,'add',[51:100]','Versicolor');
  sDiris = som_label(sDiris,'add',[101:150]','Virginica');
  
  echo on
end

%     Here are the histograms and scatter plots of the four variables.

echo off 
k=1;
for i=1:4, 
  for j=1:4, 
    if i==j, 
      subplot(4,4,k); 
      hist(sDiris.data(:,i)); title(sDiris.comp_names{i})
    elseif i<j, 
      subplot(4,4,k); 
      plot(sDiris.data(:,i),sDiris.data(:,j),'k.')
      xlabel(sDiris.comp_names{i})
      ylabel(sDiris.comp_names{j})
    end
    k=k+1;
  end
end
echo on

%     Actually, as you saw in SOM_DEMO1, most SOM Toolbox functions
%     can also handle plain data matrices, but then one is without the
%     convenience offered by component names, labels and
%     denormalization operations.


pause % Strike any key to normalize the data...





clc
%    STEP 2: DATA NORMALIZATION
%    ==========================

%    Since SOM algorithm is based on Euclidian distances, the scale of
%    the variables is very important in determining what the map will
%    be like. If the range of values of some variable is much bigger
%    than of the other variables, that variable will probably dominate
%    the map organization completely. 

%    For this reason, the components of the data set are usually
%    normalized, for example so that each component has unit
%    variance. This can be done with function SOM_NORMALIZE:

sDiris = som_normalize(sDiris,'var');

%    The function has also other normalization methods.

%    However, interpreting the values may be harder when they have
%    been normalized. Therefore, the normalization operations can be
%    reversed with function SOM_DENORMALIZE:

x = sDiris.data(1,:)

orig_x = som_denormalize(x,sDiris)

pause % Strike any key to to train the map...





clc
%    STEP 3: MAP TRAINING
%    ====================

%    The function SOM_MAKE is used to train the SOM. By default, it
%    first determines the map size, then initializes the map using
%    linear initialization, and finally uses batch algorithm to train
%    the map.  Function SOM_DEMO1 has a more detailed description of
%    the training process.

sMap = som_make(sDiris);


pause % Strike any key to continues...

%    The IRIS data set also has labels associated with the data
%    samples. Actually, the data set consists of 50 samples of three
%    species of Iris-flowers (a total of 150 samples) such that the
%    measurements are width and height of sepal and petal leaves. The
%    label associated with each sample is the species information:
%    'Setosa', 'Versicolor' or 'Virginica'.

%    Now, the map can be labelled with these labels. The best
%    matching unit of each sample is found from the map, and the
%    species label is given to the map unit. Function SOM_AUTOLABEL 
%    can be used to do this: 

sMap = som_autolabel(sMap,sDiris,'vote');

pause % Strike any key to visualize the map...





clc
%    STEP 4: VISUALIZING THE SELF-ORGANIZING MAP: SOM_SHOW
%    =====================================================

%    The basic visualization of the SOM is done with function SOM_SHOW.

colormap(1-gray)
som_show(sMap,'norm','d')

%    Notice that the names of the components are included as the
%    titles of the subplots. Notice also that the variable values
%    have been denormalized to the original range and scale.

%    The component planes ('PetalL', 'PetalW', 'SepalL' and 'SepalW')
%    show what kind of values the prototype vectors of the map units
%    have. The value is indicated with color, and the colorbar on the
%    right shows what the colors mean.

%    The 'U-matrix' shows distances between neighboring units and thus
%    visualizes the cluster structure of the map. Note that the
%    U-matrix visualization has much more hexagons that the
%    component planes. This is because distances *between* map units
%    are shown, and not only the distance values *at* the map units. 

%    High values on the U-matrix mean large distance between
%    neighboring map units, and thus indicate cluster
%    borders. Clusters are typically uniform areas of low
%    values. Refer to colorbar to see which colors mean high
%    values. In the IRIS map, there appear to be two clusters.

pause % Strike any key to continue...

%    The subplots are linked together through similar position. In
%    each axis, a particular map unit is always in the same place. For
%    example:

h=zeros(sMap.topol.msize); h(1,2) = 1;
som_show_add('hit',h(:),'markercolor','r','markersize',0.5,'subplot','all')

%    the red marker is on top of the same unit on each axis. 

pause % Strike any key to continue...



clf

clc

%    STEP 4: VISUALIZING THE SELF-ORGANIZING MAP: SOM_SHOW_ADD
%    =========================================================

%    The SOM_SHOW_ADD function can be used to add markers, labels and
%    trajectories on top of SOM_SHOW created figures. The function
%    SOM_SHOW_CLEAR can be used to clear them away.

%    Here, the U-matrix is shown on the left, and an empty grid
%    named 'Labels' is shown on the right.

som_show(sMap,'umat','all','empty','Labels')

pause % Strike any key to add labels...

%    Here, the labels added to the map with SOM_AUTOLABEL function
%    are shown on the empty grid.

som_show_add('label',sMap,'Textsize',8,'TextColor','r','Subplot',2)

pause % Strike any key to add hits...

%    An important tool in data analysis using SOM are so called hit
%    histograms. They are formed by taking a data set, finding the BMU
%    of each data sample from the map, and increasing a counter in a
%    map unit each time it is the BMU. The hit histogram shows the
%    distribution of the data set on the map.

%    Here, the hit histogram for the whole data set is calculated
%    and visualized on the U-matrix.

h = som_hits(sMap,sDiris);
som_show_add('hit',h,'MarkerColor','w','Subplot',1)

pause % Strike any key to continue...

%    Multiple hit histograms can be shown simultaniously. Here, three
%    hit histograms corresponding to the three species of Iris
%    flowers is calculated and shown. 

%    First, the old hit histogram is removed.

som_show_clear('hit',1)

%    Then, the histograms are calculated. The first 50 samples in
%    the data set are of the 'Setosa' species, the next 50 samples
%    of the 'Versicolor' species and the last 50 samples of the
%    'Virginica' species. 

h1 = som_hits(sMap,sDiris.data(1:50,:));
h2 = som_hits(sMap,sDiris.data(51:100,:));
h3 = som_hits(sMap,sDiris.data(101:150,:));

som_show_add('hit',[h1, h2, h3],'MarkerColor',[1 0 0; 0 1 0; 0 0 1],'Subplot',1)

%    Red color is for 'Setosa', green for 'Versicolor' and blue for
%    'Virginica'. One can see that the three species are pretty well
%    separated, although 'Versicolor' and 'Virginica' are slightly
%    mixed up.

pause % Strike any key to continue...



clf
clc

%    STEP 4: VISUALIZING THE SELF-ORGANIZING MAP: SOM_GRID
%    =====================================================

%    There's also another visualization function: SOM_GRID.  This
%    allows visualization of the SOM in freely specified coordinates,
%    for example the input space (of course, only upto 3D space). This
%    function has quite a lot of options, and is pretty flexible.

%    Basically, the SOM_GRID visualizes the SOM network: each unit is
%    shown with a marker and connected to its neighbors with lines.
%    The user has control over: 
%     - the coordinate of each unit (2D or 3D)
%     - the marker type, color and size of each unit
%     - the linetype, color and width of the connecting lines
%    There are also some other options.

pause % Strike any key to see some visualizations...

%    Here are four visualizations made with SOM_GRID: 
%     - The map grid in the output space.

subplot(2,2,1)
som_grid(sMap,'Linecolor','k')
view(0,-90), title('Map grid')

%     - A surface plot of distance matrix: both color and 
%       z-coordinate indicate average distance to neighboring 
%       map units. This is closely related to the U-matrix.

subplot(2,2,2)
Co=som_unit_coords(sMap); U=som_umat(sMap); U=U(1:2:size(U,1),1:2:size(U,2));
som_grid(sMap,'Coord',[Co, U(:)],'Surf',U(:),'Marker','none');
view(-80,45), axis tight, title('Distance matrix')

%     - The map grid in the output space. Three first components
%       determine the 3D-coordinates of the map unit, and the size
%       of the marker is determined by the fourth component.
%       Note that the values have been denormalized.

subplot(2,2,3)
M = som_denormalize(sMap.codebook,sMap);
som_grid(sMap,'Coord',M(:,1:3),'MarkerSize',M(:,4)*2)
view(-80,45), axis tight, title('Prototypes')

%     - Map grid as above, but the original data has been plotted
%       also: coordinates show the values of three first components
%       and color indicates the species of each sample.  Fourth
%       component is not shown.

subplot(2,2,4)
som_grid(sMap,'Coord',M(:,1:3),'MarkerSize',M(:,4)*2)
hold on
D = som_denormalize(sDiris.data,sDiris); 
plot3(D(1:50,1),D(1:50,2),D(1:50,3),'r.',...
      D(51:100,1),D(51:100,2),D(51:100,3),'g.',...
      D(101:150,1),D(101:150,2),D(101:150,3),'b.')
view(-72,64), axis tight, title('Prototypes and data')

pause % Strike any key to continue...

%    STEP 5: ANALYSIS OF RESULTS
%    ===========================

%    The purpose of this step highly depends on the purpose of the
%    whole data analysis: is it segmentation, modeling, novelty
%    detection, classification, or something else? For this reason, 
%    there is not a single general-purpose analysis function, but 
%    a number of individual functions which may, or may not, prove 
%    useful in any specific case.

%    Visualization is of course part of the analysis of
%    results. Examination of labels and hit histograms is another
%    part. Yet another is validation of the quality of the SOM (see
%    the use of SOM_QUALITY in SOM_DEMO1).

[qe,te] = som_quality(sMap,sDiris)

%    People have contributed a number of functions to the Toolbox
%    which can be used for the analysis. These include functions for 
%    vector projection, clustering, pdf-estimation, modeling,
%    classification, etc. However, ultimately the use of these
%    tools is up to you.

%    More about visualization is presented in SOM_DEMO3.
%    More about data analysis is presented in SOM_DEMO4.

echo off
warning on




