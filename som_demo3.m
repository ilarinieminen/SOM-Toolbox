
%SOM_DEMO3 Self-organizing map visualization.

% Contributed to SOM Toolbox 2.0, February 11th, 2000 by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/

% Version 1.0beta juuso 071197 
% Version 2.0beta juuso 080200 070600

clf reset;
figure(gcf)
echo on




clc
%    ==========================================================
%    SOM_DEMO3 - VISUALIZATION
%    ==========================================================

%    som_show           - Visualize map.
%    som_grid           - Visualization with free coordinates.
%
%    som_show_add       - Add markers on som_show visualization.
%    som_show_clear     - Remove markers from som_show visualization.
%    som_recolorbar     - Refresh and rescale colorbars in som_show 
%                         visualization.
%
%    som_cplane         - Visualize component/color/U-matrix plane.
%    som_pieplane       - Visualize prototype vectors as pie charts.
%    som_barplane       - Visualize prototype vectors as bar charts.
%    som_plotplane      - Visualize prototype vectors as line graphs.
%
%    pcaproj            - Projection to principal component space.
%    cca                - Projection with Curvilinear Component Analysis.
%    sammon             - Projection with Sammon's mapping.
%    som_umat           - Calculate U-matrix.
%    som_colorcode      - Color coding for the map.
%    som_normcolor      - RGB values of indexed colors.
%    som_hits           - Hit histograms for the map.

%    The basic functions for SOM visualization are SOM_SHOW and
%    SOM_GRID. The SOM_SHOW has three auxiliary functions:
%    SOM_SHOW_ADD, SOM_SHOW_CLEAR and SOM_RECOLORBAR which are used 
%    to add and remove markers and to control the colorbars.
%    SOM_SHOW actually uses SOM_CPLANE to make the visualizations.
%    Also SOM_{PIE,BAR,PLOT}PLANE can be used to visualize SOMs.

%    The other functions listed above do not themselves visualize
%    anything, but their results are used in the visualizations.
    
%    There's an important limitation that visualization functions have:
%    while the SOM Toolbox otherwise supports N-dimensional map grids, 
%    visualization only works for 1- and 2-dimensional map grids!!!

pause % Strike any key to create demo data and map...





clc
%    DEMO DATA AND MAP
%    =================

%    The data set contructed for this demo consists of random vectors
%    in three gaussian kernels the centers of which are at [0, 0, 0],
%    [3 3 3] and [9 0 0]. The map is trained using default parameters.

D1 = randn(100,3);
D2 = randn(100,3) + 3;
D3 = randn(100,3); D3(:,1) = D3(:,1) + 9;

sD = som_data_struct([D1; D2; D3],'name','Demo3 data',...
		     'comp_names',{'X-coord','Y-coord','Z-coord'});
sM = som_make(sD);

%    Since the data (and thus the prototypes of the map) are
%    3-dimensional, they can be directly plotted using PLOT3.
%    Below, the data is plotted using red 'o's and the map
%    prototype vectors with black '+'s.

plot3(sD.data(:,1),sD.data(:,2),sD.data(:,3),'ro',...
      sM.codebook(:,1),sM.codebook(:,2),sM.codebook(:,3),'k+')
rotate3d on

%    From the visualization it is pretty easy to see what the data is
%    like, and how the prototypes have been positioned. One can see
%    that there are three clusters, and that there are some prototype
%    vectors between the clusters, although there is actually no
%    data there. The map units corresponding to these prototypes
%    are called 'dead' or 'interpolative' map units.

pause % Strike any key to continue...



clc
%    VISUALIZATION OF MULTIDIMENSIONAL DATA
%    ======================================

%    Usually visualization of data sets is not this straightforward,
%    since the dimensionality is much higher than three. In principle,
%    one can embed additional information to the visualization by
%    using properties other than position, for example color, size or
%    shape.

%    Here the data set and map prototypes are plotted again, but
%    information of the cluster is shown using color: red for the
%    first cluster, green for the second and blue for the last.

plot3(sD.data(1:100,1),sD.data(1:100,2),sD.data(1:100,3),'ro',...
      sD.data(101:200,1),sD.data(101:200,2),sD.data(101:200,3),'go',...
      sD.data(201:300,1),sD.data(201:300,2),sD.data(201:300,3),'bo',...
      sM.codebook(:,1),sM.codebook(:,2),sM.codebook(:,3),'k+')
rotate3d on

%    However, this works only for relatively small dimensionality, say
%    less than 10. When the information is added this way, the
%    visualization becomes harder and harder to understand. Also, not
%    all properties are equal: the human visual system perceives
%    colors differently from position, not to mention the complex
%    rules governing perception of shape. 

pause % Strike any key to learn about linking...





clc
%    LINKING MULTIPLE VISUALIZATIONS
%    ===============================

%    The other option is to use *multiple visualizations*, so called
%    small multiples, instead of only one. The problem is then how to
%    link these visualizations together: one should be able to idetify
%    the same object from the different visualizations.

%    This could be done using, for example, color: each object has
%    the same color in each visualization. Another option is to use 
%    similar position: each object has the same position in each
%    small multiple.

%    For example, here are four subplots, one for each component and
%    one for cluster information, where color denotes the value and
%    position is used for linking. The 2D-position is derived by
%    projecting the data into the space spanned by its two greatest
%    eigenvectors.

[Pd,V,me] = pcaproj(sD.data,2);        % project the data
Pm        = pcaproj(sM.codebook,V,me); % project the prototypes
colormap(hot);                         % colormap used for values

echo off
for c=1:3, 
  subplot(2,2,c), cla, hold on
  som_grid('rect',[300 1],'coord',Pd,'Line','none',...
	   'MarkerColor',som_normcolor(sD.data(:,c)));
  som_grid(sM,'Coord',Pm,'Line','none','marker','+');
  hold off, title(sD.comp_names{c}), xlabel('PC 1'), ylabel('PC 2');
end

subplot(2,2,4), cla
plot(Pd(1:100,1),Pd(1:100,2),'ro',...
     Pd(101:200,1),Pd(101:200,2),'go',...
     Pd(201:300,1),Pd(201:300,2),'bo',...
     Pm(:,1),Pm(:,2),'k+')
title('Cluster')
echo on

pause % Strike any key to use color for linking...

%    Here is another example, where color is used for linking. On the
%    top right triangle are the scatter plots of each variable without
%    color coding, and on the bottom left triangle with the color
%    coding. In the colored figures, each data sample can be
%    identified by a unique color. Well, almost identified: there are
%    quite a lot of samples with almost the same color. Color is not as
%    precise linking method as position.

echo off 
Col = som_normcolor([1:300]',jet(300));
k=1;
for i=1:3, 
  for j=1:3, 
    if i<j, i1=i; i2=j; else i1=j; i2=i; end
    if i<j,
      subplot(3,3,k); cla
      plot(sD.data(:,i1),sD.data(:,i2),'ko')
      xlabel(sD.comp_names{i1}), ylabel(sD.comp_names{i2})
    elseif i>j,
      subplot(3,3,k); cla
      som_grid('rect',[300 1],'coord',sD.data(:,[i1 i2]),...
	       'Line','none','MarkerColor',Col);
      xlabel(sD.comp_names{i1}), ylabel(sD.comp_names{i2})
    end
    k=k+1;
  end
end
echo on

pause % Strike any key to learn about data visualization using SOM...


clc
%    DATA VISUALIZATION USING SOM
%    ============================

%    The basic visualization functions and their usage have already
%    been introduced in SOM_DEMO2. In this demo, a more structured
%    presentation is given. 

%    Data visualization techniques using the SOM can be divided to
%    three categories based on their goal:

%     1. visualization of clusters and shape of the data:
%        projections, U-matrices and other distance matrices
%
%     2. visualization of components / variables: 
%        component planes, scatter plots
%
%     3. visualization of data projections: 
%        hit histograms, response surfaces

pause % Strike any key to visualize clusters with distance matrices...



clf
clc
%    1. VISUALIZATION OF CLUSTERS: DISTANCE MATRICES
%    ===============================================

%    Distance matrices are typically used to show the cluster
%    structure of the SOM. They show distances between neighboring
%    units, and are thus closely related to single linkage clustering
%    techniques. The most widely used distance matrix technique is
%    the U-matrix. 

%    Here, the U-matrix of the map is shown (using all three
%    components in the distance calculation):

colormap(1-gray)
som_show(sM,'umat','all');

pause % Strike any key to see more examples of distance matrices...

%    The function SOM_UMAT can be used to calculate U-matrix.  The
%    resulting matrix holds distances between neighboring map units,
%    as well as the median distance from each map unit to its
%    neighbors. These median distances corresponding to each map unit
%    can be easily extracted. The result is a distance matrix using
%    median distance.

U = som_umat(sM);
Um = U(1:2:size(U,1),1:2:size(U,2));

%    A related technique is to assign colors to the map units such
%    that similar map units get similar colors.

%    Here, four clustering figures are shown: 
%     - U-matrix
%     - median distance matrix (with grayscale)
%     - median distance matrix (with map unit size)
%     - similarity coloring, made by spreading a colormap
%       on top of the principal component projection of the
%       prototype vectors

subplot(2,2,1)
h=som_cplane([sM.topol.lattice,'U'],sM.topol.msize, U(:)); 
set(h,'Edgecolor','none'); title('U-matrix')

subplot(2,2,2)
h=som_cplane(sM, Um(:));
set(h,'Edgecolor','none'); title('D-matrix (grayscale)')

subplot(2,2,3)
som_cplane(sM,'none',1-Um(:)/max(Um(:)))
title('D-matrix (marker size)')

subplot(2,2,4)
C = som_colorcode(Pm);  % Pm is the PC-projection calculated earlier
som_cplane(sM,C)
title('Similarity coloring')

pause % Strike any key to visualize shape and clusters with projections...



clf
clc
%    1. VISUALIZATION OF CLUSTERS AND SHAPE: PROJECTIONS
%    ===================================================

%    In vector projection, a set of high-dimensional data samples is
%    projected to a lower dimensional such that the distances between
%    data sample pairs are preserved as well as possible. Depending 
%    on the technique, the projection may be either linear or
%    non-linear, and it may place special emphasis on preserving
%    local distances. 

%    For example SOM is a projection technique, since the prototypes
%    have well-defined positions on the 2-dimensional map grid. SOM as
%    a projection is however a very crude one. Other projection
%    techniques include the principal component projection used
%    earlier, Sammon's mapping and Curvilinear Component Analysis
%    (to name a few). These have been implemented in functions
%    PCAPROJ, SAMMON and CCA. 

%    Projecting the map prototype vectors and joining neighboring map
%    units with lines gives the SOM its characteristic net-like look.
%    The projection figures can be linked to the map planes using
%    color coding.

%    Here is the distance matrix, color coding, a projection without
%    coloring and a projection with one. In the last projection,
%    the size of interpolating map units has been set to zero.

subplot(2,2,1)
som_cplane(sM,Um(:));
title('Distance matrix')

subplot(2,2,2)
C = som_colorcode(sM,'rgb4');
som_cplane(sM,C);
title('Color code')

subplot(2,2,3)
som_grid(sM,'Coord',Pm,'Linecolor','k');
title('PC-projection')

subplot(2,2,4)
h = som_hits(sM,sD); s=6*(h>0);
som_grid(sM,'Coord',Pm,'MarkerColor',C,'Linecolor','k','MarkerSize',s);
title('Colored PC-projection')

pause % Strike any key to visualize component planes...


clf
clc
%    2. VISUALIZATION OF COMPONENTS: COMPONENT PLANES
%    ================================================

%    The component planes visualizations shows what kind of values the
%    prototype vectors of the map units have for different vector
%    components.

%    Here is the U-matrix and the three component planes of the map.

som_show(sM)

pause % Strike any key to continue...

%    Besides SOM_SHOW and SOM_CPLANE, there are three other
%    functions specifically designed for showing the values of the 
%    component planes: SOM_PIEPLANE, SOM_BARPLANE, SOM_PLOTPLANE. 

%    SOM_PIEPLANE shows a single pie chart for each map unit.  Each
%    pie shows the relative proportion of each component of the sum of
%    all components in that map unit. The component values must be
%    positive. 

%    SOM_BARPLANE shows a barchart in each map unit. The scaling of 
%    bars can be made unit-wise or variable-wise. By default it is
%    determined variable-wise.

%    SOM_PLOTPLANE shows a linegraph in each map unit. 

M = som_normalize(sM.codebook,'range'); 

subplot(1,3,1)
som_pieplane(sM, M);
title('som\_pieplane')

subplot(1,3,2)
som_barplane(sM, M, '', 'unitwise');
title('som\_barplane')

subplot(1,3,3)
som_plotplane(sM, M, 'b');
title('som\_plotplane')

pause % Strike any key to visualize cluster properties...



clf
clc
%    2. VISUALIZATION OF COMPONENTS: CLUSTERS
%    ========================================

%    An interesting question is of course how do the values of the
%    variables relate to the clusters: what are the values of the
%    components in the clusters, and which components are the ones
%    which *make* the clusters.

som_show(sM)

%    From the U-matrix and component planes, one can easily see
%    what the typical values are in each cluster. 

pause % Strike any key to continue...

%    The significance of the components with respect to the clustering
%    is harder to visualize. One indication of importance is that on
%    the borders of the clusters, values of important variables change
%    very rabidly.

%    Here, the distance matrix is calculated with respect to each
%    variable. 

u1 = som_umat(sM,'mask',[1 0 0]'); u1=u1(1:2:size(u1,1),1:2:size(u1,2));
u2 = som_umat(sM,'mask',[0 1 0]'); u2=u2(1:2:size(u2,1),1:2:size(u2,2));
u3 = som_umat(sM,'mask',[0 0 1]'); u3=u3(1:2:size(u3,1),1:2:size(u3,2));

%    Here, the distance matrices are shown, as well as a piechart
%    indicating the relative importance of each variable in each
%    map unit. The size of piecharts has been scaled by the
%    distance matrix calculated from all components.

subplot(2,2,1)
som_cplane(sM,u1(:));
title(sM.comp_names{1})

subplot(2,2,2)
som_cplane(sM,u2(:));
title(sM.comp_names{2})

subplot(2,2,3)
som_cplane(sM,u3(:));
title(sM.comp_names{3})

subplot(2,2,4)
som_pieplane(sM, [u1(:), u2(:), u3(:)], hsv(3), Um(:)/max(Um(:)));
title('Relative importance')

%    From the last subplot, one can see that in the area where the
%    bigger cluster border is, the 'X-coord' component (red color)
%    has biggest effect, and thus is the main factor in separating
%    that cluster from the rest.

pause % Strike any key to learn about correlation hunting...


clf
clc
%    2. VISUALIZATION OF COMPONENTS: CORRELATION HUNTING
%    ===================================================

%    Finally, the component planes are often used for correlation
%    hunting. When the number of variables is high, the component
%    plane visualization offers a convenient way to visualize all
%    components at once and hunt for correlations (as opposed to
%    N*(N-1)/2 scatterplots).

%    Hunting correlations this way is not very accurate. However, it
%    is easy to select interesting combinations for further
%    investigation.

%    Here, the first and third components are shown with scatter
%    plot. As with projections, a color coding is used to link the
%    visualization to the map plane. In the color coding, size shows
%    the distance matrix information.

C = som_colorcode(sM);
subplot(1,2,1)
som_cplane(sM,C,1-Um(:)/max(Um(:)));
title('Color coding + distance matrix')

subplot(1,2,2)
som_grid(sM,'Coord',sM.codebook(:,[1 3]),'MarkerColor',C);
title('Scatter plot'); xlabel(sM.comp_names{1}); ylabel(sM.comp_names{3})
axis equal

pause % Strike any key to visualize data responses...


clf
clc
%    3. DATA ON MAP
%    ==============

%    The SOM is a map of the data manifold. An interesting question
%    then is where on the map a specific data sample is located, and
%    how accurate is that localization? One is interested in the
%    response of the map to the data sample. 

%    The simplest answer is to find the BMU of the data sample.
%    However, this gives no indication of the accuracy of the
%    match. Is the data sample close to the BMU, or is it actually
%    equally close to the neighboring map units (or even approximately
%    as close to all map units)? Sometimes accuracy doesn't really
%    matter, but if it does, it should be visualized somehow.

%    Here are different kinds of response visualizations for two
%    vectors: [0 0 0] and [99 99 99]. 
%     - BMUs indicated with labels       
%     - BMUs indicated with markers, relative quantization errors
%       (in this case, proportion between distances to BMU and 
%       Worst-MU) with vertical lines
%     - quantization error between the samples and all map units 
%     - fuzzy response (a non-linear function of quantization
%       error) of all map units

echo off
[bm,qe] = som_bmus(sM,[0 0 0; 99 99 99],'all'); % distance to all map units
[dummy,ind] = sort(bm(1,:)); d0 = qe(1,ind)'; 
[dummy,ind] = sort(bm(2,:)); d9 = qe(2,ind)'; 
bmu0 = bm(1,1); bmu9 = bm(2,1); % bmus

h0 = zeros(prod(sM.topol.msize),1); h0(bmu0) = 1; % crisp hits
h9 = zeros(prod(sM.topol.msize),1); h9(bmu9) = 1; 

lab = cell(prod(sM.topol.msize),1); 
lab{bmu0} = '[0,0,0]'; lab{bmu9} = '[99,99,99]';

hf0 = som_hits(sM,[0 0 0],'fuzzy'); % fuzzy response
hf9 = som_hits(sM,[99 99 99],'fuzzy'); 

som_show(sM,'umat',{'all','BMU'},...
	 'color',{d0,'Qerror 0'},'color',{hf0,'Fuzzy response 0'},...
	 'empty','BMU+qerror',...
	 'color',{d9,'Qerror 99'},'color',{hf9,'Fuzzy response 99'});	 
som_show_add('label',lab,'Subplot',1,'Textcolor','r');
som_show_add('hit',[h0, h9],'Subplot',4,'MarkerColor','r');
hold on
Co = som_vis_coords(sM.topol.lattice,sM.topol.msize);
plot3(Co(bmu0,[1 1]),Co(bmu0,[2 2]),[0 10*qe(1,1)/qe(1,end)],'r-')
plot3(Co(bmu9,[1 1]),Co(bmu9,[2 2]),[0 10*qe(2,1)/qe(2,end)],'r-')
view(3), axis equal
echo on

%    Here are the distances to BMU, 2-BMU and WMU:

qe(1,[1,2,end]) % [0 0 0]
qe(2,[1,2,end]) % [99 99 99]

%    One can see that for [0 0 0] the accuracy is pretty good as the
%    quantization error of the BMU is much lower than that of the
%    WMU. On the other hand [99 99 99] is very far from the map:
%    distance to BMU is almost equal to distance to WMU.

pause % Strike any key to visualize responses of multiple samples...



clc
clf
%    3. DATA ON MAP: HIT HISTOGRAMS
%    ==============================

%    One can also investigate whole data sets using the map. When the
%    BMUs of multiple data samples are aggregated, a hit histogram
%    results. Instead of BMUs, one can also aggregate for example
%    fuzzy responses.

%    The hit histograms (or aggregated responses) can then be compared
%    with each other. 

%    Here are hit histograms of three data sets: one with 50 first
%    vectors of the data set, one with 150 samples from the data
%    set, and one with 50 randomly selected samples. In the last
%    subplot, the fuzzy response of the first data set.

dlen = size(sD.data,1);
Dsample1 = sD.data(1:50,:); h1 = som_hits(sM,Dsample1); 
Dsample2 = sD.data(1:150,:); h2 = som_hits(sM,Dsample2); 
Dsample3 = sD.data(ceil(rand(50,1)*dlen),:); h3 = som_hits(sM,Dsample3); 
hf = som_hits(sM,Dsample1,'fuzzy');

som_show(sM,'umat','all','umat','all','umat','all','color',{hf,'Fuzzy'})
som_show_add('hit',h1,'Subplot',1,'Markercolor','r')
som_show_add('hit',h2,'Subplot',2,'Markercolor','r')
som_show_add('hit',h3,'Subplot',3,'Markercolor','r')

pause % Strike any key to visualize trajectories...



clc
clf
%    3. DATA ON MAP: TRAJECTORIES
%    ============================

%    A special data mapping technique is trajectory. If the samples
%    are ordered, forming a time-series for example, their response on
%    the map can be tracked. The function SOM_SHOW_ADD can be used to
%    show the trajectories in two different modes: 'traj' and 'comet'.

%    Here, a series of data points is formed which go from [8,0,0]
%    to [2,2,2]. The trajectory is plotted using the two modes.

Dtraj = [linspace(9,2,20); linspace(0,2,20); linspace(0,2,20)]';
T = som_bmus(sM,Dtraj);

som_show(sM,'comp',[1 1]);
som_show_add('traj',T,'Markercolor','r','TrajColor','r','subplot',1);
som_show_add('comet',T,'MarkerColor','r','subplot',2);

%    There's also a function SOM_TRAJECTORY which lauches a GUI
%    specifically designed for displaying trajectories (in 'comet'
%    mode).

pause % Strike any key to learn about color handling...




clc
clf
%    COLOR HANDLING
%    ==============

%    Matlab offers flexibility in the colormaps. Using the COLORMAP
%    function, the colormap may be changed. There are several useful
%    colormaps readily available, for example 'hot' and 'jet'. The
%    default number of colors in the colormaps is 64. However, it is
%    often advantageous to use less colors in the colormap. This way
%    the components planes visualization become easier to interpret.

%    Here the three component planes are visualized using the 'hot'
%    colormap and only three colors.

som_show(sM,'comp',[1 2 3])
colormap(hot(3));
som_recolorbar 

pause % Press any key to change the colorbar labels...

%     The function SOM_RECOLORBAR can be used to reconfigure
%     the labels beside the colorbar. 

%     Here the colorbar of the first subplot is labeled using labels
%     'small', 'medium' and 'big' at values 0, 1 and 2.  For the
%     colorbar of the second subplot, values are calculated for the
%     borders between colors.

som_recolorbar(1,{[0 4 9]},'',{{'small','medium','big'}});
som_recolorbar(2,'border','');

pause % Press any key to learn about SOM_NORMCOLOR...

%     Some SOM Toolbox functions do not use indexed colors if the
%     underlying Matlab function (e.g. PLOT) do not use indexed
%     colors. SOM_NORMCOLOR is a convenient function to simulate
%     indexed colors: it calculates fixed RGB colors that
%     are similar to indexed colors with the specified colormap. 

%     Here, two SOM_GRID visualizations are created. One uses the
%     'surf' mode to show the component colors in indexed color
%     mode, and the other uses SOM_NORMALIZE to do the same. 

clf
colormap(jet(64))
subplot(1,2,1)
som_grid(sM,'Surf',sM.codebook(:,3));
title('Surf mode')

subplot(1,2,2)
som_grid(sM,'Markercolor',som_normcolor(sM.codebook(:,3)));
title('som\_normcolor')

pause % Press any key to visualize different map shapes...



clc
clf
%    DIFFERENT MAP SHAPES
%    ====================

%    There's no direct way to visualize cylinder or toroid maps. When
%    visualized, they are treated exactly as if they were sheet
%    shaped. However, if function SOM_UNIT_COORDS is used to provide
%    unit coordinates, then SOM_GRID can be used to visualize these
%    alternative map shapes.

%    Here the grids of the three possible map shapes (sheet, cylinder
%    and toroid) are visualized. The last subplot shows a component 
%    plane visualization of the toroid map.

Cor = som_unit_coords(sM.topol.msize,'hexa','sheet');
Coc = som_unit_coords(sM.topol.msize,'hexa','cyl');
Cot = som_unit_coords(sM.topol.msize,'hexa','toroid');

subplot(2,2,1)
som_grid(sM,'Coord',Cor,'Markersize',3,'Linecolor','k');
title('sheet'), view(0,-90), axis tight, axis equal

subplot(2,2,2)
som_grid(sM,'Coord',Coc,'Markersize',3,'Linecolor','k');
title('cylinder'), view(5,1), axis tight, axis equal

subplot(2,2,3)
som_grid(sM,'Coord',Cot,'Markersize',3,'Linecolor','k');
title('toroid'), view(-100,0), axis tight, axis equal

subplot(2,2,4)
som_grid(sM,'Coord',Cot,'Surf',sM.codebook(:,3));
colormap(jet), colorbar
title('toroid'), view(-100,0), axis tight, axis equal

echo off
