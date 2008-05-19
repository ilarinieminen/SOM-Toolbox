function vis_trajgui(trajStruct,arg)

% VIS_TRAJGUI subfuntion for SOM_TRAJECTORY
%
% This function is the actual GUI called by SOM_TRAJECTORY
% function. 
%
% See also SOM_TRAJECTORY.

% Contributed code to SOM Toolbox 2.0, February 11th, 2000 by Juha Parhankangas
% Copyright (c) by Juha Parhankangas.
% http://www.cis.hut.fi/projects/somtoolbox/        

% Version 2.0beta juha 180699

if nargin == 1

  sM_h=trajStruct.figure;

  if size(trajStruct.bmus,1) ~= 1 & size(trajStruct.bmus,2) ~= 1
    fuzzy_traj(trajStruct,[]);
    return;
  end
  
  
  if size(trajStruct.bmus,1) == 1 | size(trajStruct.bmus,2) == 1

    udata.bmus = trajStruct.bmus;
    udata.a_h=[findobj(get(sM_h,'Children'),'Tag','Uplane');...
	       findobj(get(sM_h,'Children'),'Tag','Cplane')];
    udata.sM_h=trajStruct.figure;
    udata.traj=[];
    data1 = trajStruct.primary_data;
    if ~isempty(trajStruct.primary_names)
      names=trajStruct.primary_names;
    else
      for i=1:size(data1,2)
	names{i,1}=sprintf('Var%d',i);
      end
    end

    udata.lattice=trajStruct.lattice;
    form = 0.7*vis_patch(udata.lattice);
    udata.msize = trajStruct.msize;


    %%%%%%%%%%%%%%%%%%%%%%%%
    %
    % forming a patch object, which is placed above every component plane
    %
    
    
    l = size(form,1);
    
    nx = repmat(form(:,1),1,prod(udata.msize));
    ny = repmat(form(:,2),1,prod(udata.msize));
    
    x=reshape(repmat(1:udata.msize(2),l*udata.msize(1),1),l,prod(udata.msize));
    y=repmat(repmat(1:udata.msize(1),l,1),1,udata.msize(2));
    
    if strcmp(udata.lattice,'hexa')
      t = find(~rem(y(1,:),2));
      x(:,t)=x(:,t)+.5;
    end
    x=x+nx;
    y=y+ny;
    
    colors=reshape(ones(prod(udata.msize),1)*[NaN NaN NaN],...
		   [1 prod(udata.msize) 3]);
    
    set(0,'CurrentFigure',udata.sM_h);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %
    % drawing patch
    %
    % caxis -commands keep the colormap of the original patch unchanged.
    %
    
    for i=1:length(udata.a_h)
      udata.real_patch(i)=get(udata.a_h(i),'Children');
      set(udata.real_patch(i),'ButtonDownFcn',...
			'vis_trajgui([],''click'')');
      subplot(udata.a_h(i));
      v=caxis;
      udata.tmp_patch(i)=patch(x,y,colors,'EdgeColor','none',...
			       'ButtonDownFcn',...
			       'vis_trajgui([],''click'')',...
			       'Tag','TmpPatch');
      caxis(v);
    end
    
    %%%%%%%%%%%%%%%%%%%%
    
    
    
    udata.length_of_traj=length(trajStruct.size);
    udata.size=trajStruct.size;
    udata.color=trajStruct.color;
    udata.poly.x=[];
    udata.poly.y=[];
    udata.poly.h=[];
    udata.new_marks=[];
    udata.all_marks=[];
    udata.d_mark2=[];
    udata.fig1 = figure;
    set(udata.fig1,'KeyPressFcn','vis_trajgui([],''key'')',...
		   'Name','Primary Data');
    
    
    %%%%%%%%%%%%%%%%%%%%
    %
    % making the 'Tools' -menu
    %
    
    udata.m_i=uimenu(udata.fig1,'Label','Trajectoy T&ools');
    udata.m_i(2)=uimenu(udata.m_i,'Label','&Remove Trajectory',...
			'Callback',...
			'vis_trajgui([],''remove_traj'')');
    udata.m_i(3)=uimenu(udata.m_i(1),'Label','&Dye Nodes',...
			'Callback',...
			'vis_trajgui([],''dye_gui'')');
    udata.m_i(4)=uimenu(udata.m_i(1),'Label','&Clear Markers',...
			'Callback',...
			'vis_trajgui([],''clear'')');
    udata.m_i(5)=uimenu(udata.m_i(1),'Label','&Save',...
			'Callback',...
			'vis_trajgui([],''save'')');
    udata.m_i(6)=uimenu(udata.m_i(1),'Label','&Load',...
			'Callback',...
			'vis_trajgui([],''load'')');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % drawing data components to the figure ....
    %
    % 
    
    if nargin < 5 | isempty(comps) | (isstr(comps) & strcmp(comps,'all'))
      comps = 1:size(data1,2);
    end
    
    x=1:size(data1,1);
    
    for i=1:length(comps)
      subplot(length(comps),1,i);
      udata.h(i)=gca;
      udata.d_mark(i).h=[];
      
      udata.d(i)=plot(x,data1(:,comps(i)),...
		      'ButtonDownFcn',...
		      'vis_trajgui([],''line_down'')');        
      set(gca,'XLim',[1 size(data1,1)],...
	      'XTick',[],...	      
	      'ButtonDownFcn','vis_trajgui([],''line_down'')'); %,...
	      %'YLim',[min(data1(:,comps(i))) max(data1(:,comps(i)))]);
	      
      ylabel(names{comps(i)});
      hold on;
      ymin=get(udata.h(i),'YLim');
      pos=mean(get(udata.h(i),'XLim'));
      udata.l(i) = line([pos pos],[ymin(1) ymin(2)],...
			'Color','red',...
			'ButtonDownFcn',...
			'vis_trajgui([],''down'')');
    end
    udata.text1=[];
    
    udata.fig2=[];
    udata.h2=[];
    udata.l2=[];
    udata.text2=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % ... and to the figure 2.
    %
    
    
    if ~isempty(trajStruct.secondary_data)
      data2=trajStruct.secondary_data;
      if isempty(trajStruct.secondary_names)
	for i=1:size(data1,2)
	  names2{i,1}=sprintf('Var%d',i);
	end
      else
	names2=trajStruct.secondary_names;
      end
      
      udata.fig2 = figure;
      set(udata.fig2,'Name','Secondary Data');
      set(udata.fig2,'KeyPressFcn',...
		     'vis_trajgui([],''key'')');
      for i=1:size(data2,2)
	subplot(size(data2,2),1,i);
	udata.h2(i) = gca;
	udata.d_mark2(i).h=[];
	udata.d2(i) = plot(x,data2(:,i),...
			   'ButtonDownFcn',...
			   'vis_trajgui([],''line_down'')');
	set(gca,'XLim',[1 size(data1,1)],'XTick',[],'ButtonDownFcn',...
		'vis_trajgui([],[],[],[],[],[],''line_down'')');
	ylabel(names2{i});
	hold on;
	ymin = get(udata.h2(i),'YLim');
	pos = mean(get(udata.h2(i),'XLim'));
	udata.l2(i) = line([pos pos],ymin,'Color','red',...
			   'ButtonDownFcn','vis_trajgui([],''down'')');
      end
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    set(udata.fig1,'UserData',udata);
    if ~isempty(udata.fig2);
      tmp.fig1=udata.fig1;
      set(udata.fig2,'UserData',tmp);
    end
    tmp=get(udata.sM_h,'UserData');
    tmp.fig1=udata.fig1;
    set(udata.sM_h,'UserData',tmp);  
    set_numbers(round(pos)); 
    return;
  end
  
end

%%%%%%%%%%%%%%%%%
%
% if figures have been drawn, the only function calls that may exist
% are the ones that change the state of the application.
%
udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

switch arg
 case 'fuzzy'
  fuzzy_traj(tS,[]);
  return;
 case 'move_fuzzy'
  fuzzy_traj([],'move');
  return;
 case 'remove_traj'
  remove_traj;
  return;
 case 'line_down'
  line_bdf('down');
  return;
 case 'line_drag'
  line_bdf('drag');
  return;
 case 'line_up'
  line_bdf('up');
  return;
 case 'dye_gui';
  color_gui(udata.fig1);
  return;
 case {'dye','cyan','magenta','yellow','red','green','blue','white','grey'}
  dye_nodes(arg);
  return;
 case 'clear'
  clear_markers;
  return;
 case 'key'
  key_bdf;
  return;   
 case 'click'
  click;
  return;
 case 'save'
  save_data;
  return;
 case 'load'
  load_data;
  return;
end


%%%%%%%%%%
%
% lines in the data figure(s) are dragged ...
%

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');


lims=get(gca,'XLim');
x = getfield(get(gca,'CurrentPoint'),{1});  % the location of the line

if x < lims(1)
  x=lims(1);
elseif x > lims(2)
  x=lims(2);
end

old = gcf;

switch arg
 case 'down',...
      
  % mouse button is pressed down above the line
  
  set(gcf,'WindowButtonMotionFcn','vis_trajgui([],''drag'')');
  set(gcf,'WindowButtonUpFcn','vis_trajgui([],''up'')');
  set(udata.l,'EraseMode','xor');
  delete(udata.text1);
  if ~isempty(udata.l2);
    set(udata.l2,'EraseMode','xor');
    delete(udata.text2);
  end
  set(gcf,'Pointer','crosshair');
  
 case 'drag'
  % change the location of the lines
  
  set(0,'CurrentFigure',udata.fig1);
  set(udata.l,'XData',[x x]);
  if ~isempty(udata.fig2)
    set(0,'CurrentFigure',udata.fig2);
    set(udata.l2,'XData',[x x]);
  end
  draw_traj(round(x));
  set(0,'CurrentFigure',old);
 case 'up'
  
  % draw trajectory and set figure to the normal state.
  
  set(udata.l,'EraseMode','normal');
  set(gcf,'Pointer','arrow','WindowButtonMotionFcn','',...
	  'WindowButtonUpFcn','');
  draw_traj(round(x));
  set_numbers(round(x));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_traj(point)

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');
color=udata.color;
eMode='normal';
if isstr(udata.color) & strcmp(udata.color,'xor')
  eMode='xor';
  color='black';
end
ind=udata.bmus(point);
[i j] = ind2sub(udata.msize,ind);
if ~mod(i,2)
  j=j+0.5;
end
old = gcf;
set(0,'CurrentFigure',udata.sM_h);
hold on;
if isempty(udata.traj) | length(udata.traj.h) ~= length(udata.a_h)
  
  % trajectory does not exist
  
  for i=1:length(udata.a_h)
    subplot(udata.a_h(i));
    hold on;
    new.h = plot(j,i,'Color',color,'EraseMode',eMode,'LineStyle','none');
    udata.traj.h(i)=new;
    udata.traj.j=j;
    udata.traj.i=i;
  end
else
  if length(udata.traj.j) == udata.length_of_traj
    % if the length of trajectory == ..., 
    udata.traj.j(1) = [];         % the first (the oldest) coordinate pair
    udata.traj.i(1) = [];         % is removed.
  end
  udata.traj.j=[udata.traj.j;j]; % the new point is added to the
  udata.traj.i=[udata.traj.i;i]; % end of coordinate vectors (i and j)
  for i=1:length(udata.a_h)
    subplot(udata.a_h(i));            % remove the existing trajectory
    delete(udata.traj.h(i).h);          % and plot the new one.
    for j=1:length(udata.traj.j)
      udata.traj.h(i).h(j)=plot(udata.traj.j(j),udata.traj.i(j),...
				'Color',color,...
				'EraseMode',eMode,'Marker','o','LineWidth',2,...
				'MarkerSize',udata.size(udata.length_of_traj-j+1),...
				'LineStyle','none');
    end
  end
end 
set(0,'CurrentFigure',udata.fig1);
set(udata.fig1,'UserData',udata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_numbers(x);

% This function writes the numbers beside of the pointer lines


udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
p = ylim(1) + 0.9*(ylim(2)-ylim(1));

old = gcf;
set(0,'CurrentFigure',udata.fig1);


for i=1:length(udata.h)
  subplot(udata.h(i));
  
  % check if the text is placed to the left side of the line...
  
  if abs(x-xlim(1)) > (abs(x-xlim(2)))
    udata.text1(i)=text(x-1,p,sprintf('%d ->',x),...
                        'VerticalAlignment','top',...
                        'HorizontalAlignment','right',...
                        'FontWeight','demi');
  else
    
    %  or to the right side.
    
    udata.text1(i)=text(x+1,p,sprintf('<- %d',x),...
			'VerticalAlignment','top',...
			'FontWeight','demi');
  end
end

if ~isempty(udata.fig2)
  set(0,'CurrentFigure',udata.fig2);
  
  for i=1:length(udata.h2)
    subplot(udata.h2(i));
    
    if abs(x-xlim(1)) > (abs(x-xlim(2)))
      udata.text2(i)=text(x-1,p,sprintf('%d ->',x),...
			  'VerticalAlignment','top',...
			  'HorizontalAlignment','right',...
			  'FontWeight','demi');
    else
      udata.text2(i)=text(x+1,p,sprintf('<- %d',x),...
			  'VerticalAlignment','top',...
			  'FontWeight','demi');
    end
  end
end

set(0,'CurrentFigure',old);
set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function remove_traj()

% delete trajectory -object from every component plane.

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

if isempty(udata.traj)
  return;
end


for i=1:length(udata.traj.h)
  delete(udata.traj.h(i).h);
end

udata.traj=[];
set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function line_bdf(arg)

% this function takes care of action when region is selected in the
% data figure.


udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');
xlim=get(gca,'XLim');

if ~(any(strcmp('THIS',fieldnames(udata))))
  p = getfield(get(gca,'CurrentPoint'),{1});
else
  p = getfield(get(udata.THIS,'CurrentPoint'),{1});
end

if p < xlim(1)
  p = xlim(1);
elseif p > xlim(2)
  p = xlim(2);
end



switch arg
 case 'down'
  
  % the mouse button is pressed down, the pointer lines are drawn.
  % and the state of the figure is changed.
  
  udata.THIS=gca;
  set(gcf,'WindowButtonMotionFcn',...
	  'vis_trajgui([],''line_drag'')',...
	  'WindowButtonUpFcn','vis_trajgui([],''line_up'')');
  udata.start_p=p;
  
  old = gcf;
  set(0,'CurrentFigure',udata.fig1);
  
  for i=1:length(udata.h)
    subplot(udata.h(i));
    udata.t_line.h(i)=line([p p],get(gca,'YLim'),'Color','red');
    udata.t_line.h2(i)=line([p p],get(gca,'YLim'),'Color','red',...
			    'EraseMode','xor');
  end
  if ~isempty(udata.h2)
    set(0,'CurrentFigure',udata.fig2);
    for i=1:length(udata.h2)
      subplot(udata.h2(i));
      udata.t_line2.h(i)=line([p p],get(gca,'YLim'),'Color','red');
      udata.t_line2.h2(i)=line([p p],get(gca,'YLim'),'Color','red',...
			       'EraseMode','xor');
    end
  end
  
 case 'drag'
  
  % change the position of the pointer lines
  
  old = gcf;
  set(0,'CurrentFigure',udata.fig1);
  set(udata.t_line.h2,'XData',[p p]);
  if ~isempty(udata.fig2)
    set(0,'CurrentFigure',udata.fig2);
    set(udata.t_line2.h2,'XData',[p p]);
  end
  set(0,'CurrentFigure',old);
 case 'up'
  
  
  % sort the 'points' -vector and draw the markers to the data and nodes
  
  points=sort([round(udata.start_p) round(p)]);
  draw_markers(points(1):points(2));
  udata=get(udata.fig1,'UserData');
  udata.new_marks=unique([udata.new_marks ;(points(1):points(2))']);
  delete([udata.t_line.h2 udata.t_line.h]);
  if ~isempty(udata.fig2)
    delete([udata.t_line2.h2 udata.t_line2.h]);
  end
  set(get(udata.THIS,'Parent'),'WindowButtonMotionFcn','',...
		    'WindowButtonUpFcn','');
  udata=rmfield(udata,'THIS');
end

set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_markers(x);

plot2data(x);
plot2plane(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot2data(x);

% plot black markers to the data figure(s)

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

old = gcf;

set(0,'CurrentFigure',udata.fig1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if there already exist points in the positions that are members
% of the set x, then the old points are removed from data figures...

for i=1:length(udata.d_mark(1).h)
  tmp1 = get(udata.d_mark(1).h(i),'XData');
  tmp2 = setdiff(tmp1,x);
  if length(tmp1) ~= length(tmp2)
    inds=[];
    for j=1:length(tmp2);
      inds=[inds find(tmp2(j)==tmp1)];
    end
    for j=1:length(udata.d_mark)
      ydata=getfield(get(udata.d_mark(j).h(i),'YData'),{inds});
      set(udata.d_mark(j).h(i),'XData',tmp2,'YData',ydata);
    end
    if ~isempty(udata.fig2)
      for j=1:length(udata.d_mark2)
        ydata=getfield(get(udata.d_mark2(j).h(i),'YData'),{inds});
        set(udata.d_mark2(j).h(i),'XData',tmp2,'YData',ydata);
      end
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ... and the new ones are plotted.

for i=1:length(udata.h)
  subplot(udata.h(i));
  h=plot(x,getfield(get(udata.d(i),'YData'),{x}),'oblack',...
	 'ButtonDownFcn',...
	 'vis_trajgui([],''line_down'')');
  udata.d_mark(i).h=[udata.d_mark(i).h;h];
end

if ~isempty(udata.h2)
  set(0,'CurrentFigure',udata.fig2);
  
  for i=1:length(udata.h2)
    subplot(udata.h2(i));
    h=plot(x,getfield(get(udata.d2(i),'YData'),{x}),'oblack',...
	   'ButtonDownFcn',...
	   'vis_trajgui([],''line_down'')');
    udata.d_mark2(i).h=[udata.d_mark2(i).h;h];
  end
end

set(0,'CurrentFigure',old);
set(udata.fig1,'UserData',udata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot2plane(x);

% sets markers to the component planes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% actually new markers are never plotted, but the color of the patch
% lying above the original component plane patch is changed black in
% the right positions.

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

udata.new_marks=unique([udata.new_marks ;x']);

for i=1:length(udata.a_h)
  col=get(udata.tmp_patch(i),'FaceVertexCData');
  if length(size(col)) == 3
    col = reshape(col,[size(col,1) 3]);
  end
  for j=1:length(udata.new_marks)
    col(udata.bmus(udata.new_marks(j)),:)=[0 0 0];
  end
  set(udata.tmp_patch(i),'FaceVertexCData',col);
end

set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function color_gui(fig1)

% construct the graphical user interface for changing the color of the
% black (marked) nodes.


udata=get(fig1,'UserData');

a = figure('Color',[0.8 0.8 0.8], ...
	   'Name','Colors', ...
	   'PaperType','a4letter', ...
	   'Position',[518 456 120 311], ...
	   'Tag','Fig1');

udata.c_struct.fig=a;

b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.701961 0.701961 0.701961], ...
	      'Position',[0.0700415 0.28956 0.830492 0.594566], ...
	      'Style','frame', ...
	      'Tag','Frame1');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'Position',[0.100059 0.301143 0.770456 0.571399], ...
	      'Style','frame', ...
	      'Tag','Frame2');

b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''cyan'')', ...
	      'Position',[0.130077 0.795326 0.170101 0.0617729], ...
	      'Style','radiobutton', ...
	      'Tag','cyan', ...
	      'Value',1);


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''magenta'')', ...
	      'Position',[0.130077 0.733553 0.170101 0.057912], ...
	      'Style','radiobutton', ...
	      'Tag','magenta');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''yellow'')', ...
	      'Position',[0.130077 0.664059 0.170101 0.0617729], ...
	      'Style','radiobutton', ...
	      'Tag','yellow');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''red'')', ...
	      'Position',[0.130077 0.590703 0.170101 0.0617729], ...
	      'Style','radiobutton', ...
	      'Tag','red');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''green'')', ...
	      'Position',[0.130077 0.525068 0.170101 0.057912], ...
	      'Style','radiobutton', ...
	      'Tag','green');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''blue'')', ...
	      'Position',[0.130077 0.455575 0.170101 0.0617729], ...
	      'Style','radiobutton', ...
	      'Tag','blue');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''white'')', ...
	      'Position',[0.130077 0.38608 0.170101 0.0617729], ...
	      'Style','radiobutton', ...
	      'Tag','white');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''grey'')', ...
	      'Position',[0.130077 0.320447 0.170101 0.057912], ...
	      'Style','radiobutton', ...
	      'Tag','grey');


b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.795326 0.470278 0.0501905], ...
	      'String','Cyan', ...
	      'Style','text', ...
	      'Tag','StaticText1');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.733553 0.520308 0.0463296], ...
	      'String','Magenta', ...
	      'Style','text', ...
	      'Tag','StaticText2');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.664059 0.470278 0.0501905], ...
	      'String','Yellow', ...
	      'Style','text', ...
	      'Tag','StaticText3');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.590703 0.470278 0.0501905], ...
	      'String','Red', ...
	      'Style','text', ...
	      'Tag','StaticText4');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.525068 0.470278 0.0463296], ...
	      'String','Green', ...
	      'Style','text', ...
	      'Tag','StaticText5');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.455575 0.470278 0.0463296], ...
	      'String','Blue', ...
	      'Style','text', ...
	      'Tag','StaticText6');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.38608 0.470278 0.0501905], ...
	      'String','White', ...
	      'Style','text', ...
	      'Tag','StaticText7');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.32019 0.320447 0.470278 0.0463296], ...
	      'String','Grey', ...
	      'Style','text', ...
	      'Tag','StaticText8');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Position',[0.0700415 0.146711 0.830492 0.135128], ...
	      'Style','frame', ...
	      'Tag','Frame3');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'Position',[0.100059 0.158293 0.770456 0.111963], ...
	      'Style','frame', ...
	      'Tag','Frame4');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[0.8 0.8 0.8], ...
	      'FontWeight','demi', ...
	      'HorizontalAlignment','left', ...
	      'Position',[0.130077 0.177597 0.270833 0.0617729], ...
	      'String','RGB', ...
	      'Style','text', ...
	      'Tag','StaticText9');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'BackgroundColor',[1 1 1], ...
	      'Position',[0.410243 0.173736 0.420249 0.0810768], ...
	      'Style','edit', ...
	      'Tag','EditText1');

udata.c_struct.RGB=b;

b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','vis_trajgui([],''dye'')', ...
	      'FontWeight','demi', ...
	      'Position',[0.0700415 0.0270256 0.360214 0.0772162], ...
	      'String','OK', ...
	      'Tag','Pushbutton1');
b = uicontrol('Parent',a, ...
	      'Units','normalized', ...
	      'Callback','close gcf', ...
	      'FontWeight','demi', ...
	      'Position',[0.54032 0.0270256 0.360214 0.0772162], ...
	      'String','Close', ...
	      'Tag','Pushbutton2');

udata.c_struct.color=[0 1 1];

tmp.fig1=fig1;
set(a,'UserData',tmp);
set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dye_nodes(arg)

% takes care of the action, when radiobuttons are pressed 
% (or the RGB value is set) in the color_gui -figure.
% It also handles the starting of dying nodes and plots.

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');


switch arg
 case {'cyan','magenta','yellow','red','green','blue','white','grey'}
  h=findobj(get(gcf,'Children'),'Style','radiobutton');
  set(h,'Value',0);
  set(gcbo,'Value',1);
end


switch arg
 case 'cyan'
  RGB = [0 1 1];
 case 'magenta'
  RGB = [1 0 1];
 case 'yellow'
  RGB = [1 1 0];
 case 'red'
  RGB = [1 0 0];
 case 'green'
  RGB = [0 1 0];
 case 'blue'
  RGB = [0 0 1];
 case 'white'
  RGB = [1 1 1];
 case 'grey'
  RGB = [0.4 0.4 0.4];
 case 'dye'
  
  RGB = get(udata.c_struct.RGB,'String');
  if isempty(RGB)
    dye;
    return;
  else
    str1='The value of RGB must be vector containing three scalars';
    str2='between 0 and 1.';
    color = str2num(RGB);
    set(udata.c_struct.RGB,'String','');
    if isempty(color)
      close gcf;
      udata=rmfield(udata,'c_struct');
      set(udata.fig1,'UserData',udata);
      errordlg([{str1};{str2}]);
      return;
    end
    if ~all([1 3] == size(color)) & ~all([3 1] == size(color))
      close gcf;
      errordlg([{str1};{str2}]);
      udata=rmfield(udata,'c_struct',udata);
      set(udata.fig1,'UserData',udata);
      return;
    end
    if ~isempty(cat(2,find(color>1),find(color<0)))
      close gcf
      errordlg([{str1};{str2}]);
      udata=rmfield(udata,'c_struct',udata);
      set(udata.fig1,'UserData',udata);
      return;
    end
    udata.c_struct.color=color;
    set(udata.fig1,'UserData',udata);
    dye;
    return;
  end
end

udata.c_struct.color=RGB;
set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dye()

% dyes black markers in the component planes and in the data figures

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

inds=unique([udata.all_marks ; udata.new_marks]);


for i=1:length(udata.d_mark);
  for j=1:length(udata.d_mark(i).h)
    if all(get(udata.d_mark(i).h(j),'Color') == [0 0 0])
      set(udata.d_mark(i).h(j),'Color',udata.c_struct.color);
    end
  end
end

if ~isempty(udata.fig2);
  for i=1:length(udata.d_mark2)
    for j=1:length(udata.d_mark2(i).h)
      if all(get(udata.d_mark2(i).h(j),'Color') == [0 0 0])
        set(udata.d_mark2(i).h(j),'Color',udata.c_struct.color);
      end
    end
  end
end


for i=1:length(udata.a_h)
  col=get(udata.tmp_patch(i),'FaceVertexCData');
  for j=1:length(udata.new_marks)
    col(udata.bmus(udata.new_marks(j)),:)=udata.c_struct.color;
  end
  set(udata.tmp_patch(i),'FaceVertexCData',col);
end


udata.all_marks=unique([udata.all_marks;udata.new_marks]);
udata.new_marks=[];
close gcf;
udata=rmfield(udata,'c_struct');
set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clear_markers()

% removes markers from the componentplanes and the data figure(s).

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

for i=1:length(udata.d_mark)
  delete(udata.d_mark(i).h);
  udata.d_mark(i).h=[];
end

for i=1:length(udata.d_mark2)
  delete(udata.d_mark2(i).h);
  udata.d_mark2(i).h=[];
end

col=NaN*get(udata.tmp_patch(1),'FaceVertexCData');
col=reshape(col,[size(col,1) 3]);

for i=1:length(udata.tmp_patch)
  set(udata.tmp_patch(i),'FaceVertexCData',col);
end

udata.new_marks=[];
udata.all_marks=[];


if any(strcmp('c_struct',fieldnames(udata)))
  udata=rmfield(udata,'c_struct');
end

set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function key_bdf

% moves trajectory and pointer lines, when either of 
% the keys '>' or '<' is pressed.

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

key=get(gcbo,'CurrentCharacter');

% The easiest way to get a new coordinates is to get them from the texts...
% The texts are either '<- x' or 'x ->' 

x=get(udata.text1(1),'String');
x=str2num(x(4:length(x)));

if isempty(x)
  x=get(udata.text1(1),'String');
  x=str2num(x(1:length(x)-3));
end

switch(key)
 case '<'
  if x ~= 1
    x= x-1;
  end
 case '>'
  if x ~= getfield(get(get(udata.text1(1),'Parent'),'XLim'),{2}) 
    x = x+1;
  end
 otherwise
  return;
end

set(udata.l,'XData',[x x]);
if ~isempty(udata.fig2)
  set(udata.l2,'XData',[x x]);
end

delete(udata.text1);
delete(udata.text2);

set_numbers(x);
draw_traj(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function click()


switch get(gcf,'SelectionType')
 case 'open'
  return;
 case {'normal','alt'}
  draw_poly;
 case 'extend'
  click_node;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function click_node()

% takes care of the action, when the middle mouse button is
% pressed (mouse pointer is above some component plane).

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');
new_marks=[];

old=gcf;

NEW=0;
AGAIN = 0;


coords=get(gca,'CurrentPoint');
row=round(coords(1,2));

if strcmp(udata.lattice,'hexa') & ~mod(row,2)
  col = round(coords(1,1) - 0.5);
else
  col = round(coords(1,1));
end

ind = sub2ind(udata.msize,row,col);
new_marks=find(udata.bmus==ind);

if strcmp(get(gcbo,'Tag'),'TmpPatch');
  
  % if the callback is made via temporary patch object, node is marked
  % (node is black) => the mark is to be removed 
  
  node_color = getfield(get(gcbo,'FaceVertexCData'),{ind,[1:3]});
  AGAIN = 1;
end



for i=1:length(udata.tmp_patch)
  color = get(udata.tmp_patch(i),'FaceVertexCData');
  if length(size(color)) ~= 2
    color = reshape(color,[size(color,1) 3]);
  end
  if all(isnan(color(ind,:)))
    NEW=1;
    color(ind,:)=[0 0 0];
  else
    color(ind,:)=[NaN NaN NaN];
  end
  set(udata.tmp_patch(i),'FaceVertexCData',color);
end

set(0,'CurrentFigure',udata.fig1);

for j=1:length(udata.h)
  subplot(udata.h(j));
  if NEW
    y=getfield(get(udata.d(j),'YData'),{new_marks});
    udata.d_mark(j).h=[udata.d_mark(j).h;plot(new_marks,y,'Color',[0 0 0],...
                                              'LineStyle','none',...
                                              'Marker','o')];
  end
end


if ~isempty(udata.fig2)
  set(0,'CurrentFigure',udata.fig2);
  for j=1:length(udata.h2);
    subplot(udata.h2(j));
    if NEW
      y=getfield(get(udata.d2(j),'YData'),{new_marks});
      udata.d_mark2(j).h=[udata.d_mark2(j).h;plot(new_marks,y,...
						  'LineStyle','none',...
						  'Color','black',...
						  'Marker','o')];
    end
  end
end

if NEW
  udata.new_marks=[udata.new_marks; new_marks];
end

if AGAIN
  
  % find marks from the data that map to the clicked node. if the color
  % of the mark(s) is the same as the node's color, remove mark(s), else
  % let mark be unchanged.
  
  for i=1:length(udata.d_mark(1).h)
    if all(node_color==get(udata.d_mark(1).h(i),'Color'))
      tmp1 = get(udata.d_mark(1).h(i),'XData');
      tmp2 = setdiff(tmp1,new_marks);
      if length(tmp1) ~= length(tmp2)
        inds=[];
        for j=1:length(tmp2);
          inds=[inds find(tmp2(j)==tmp1)];
        end
        for j=1:length(udata.d_mark)
          ydata=getfield(get(udata.d_mark(j).h(i),'YData'),{inds});
          set(udata.d_mark(j).h(i),'XData',tmp2,'YData',ydata);
        end
        if ~isempty(udata.fig2)
          for j=1:length(udata.d_mark2)
            ydata=getfield(get(udata.d_mark2(j).h(i),'YData'),{inds});
            set(udata.d_mark2(j).h(i),'XData',tmp2,'YData',ydata);
          end
        end
      end  
    end
  end
  udata.new_marks=setdiff(udata.new_marks, new_marks);
  udata.all_marks=setdiff(udata.all_marks,new_marks);
end

set(udata.fig1,'UserData',udata);
set(0,'CurrentFigure',old);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_poly()

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

if isempty(udata.poly.x)
  if strcmp(get(gcf,'SelectionType'),'alt')
    return;
  end
  udata.poly.THIS = gca;
end

% 'THIS' indicates what was the axes where the polygon was meant to
% drawn. It is not possible to add points, that lie in another axes, to the
% polygon.


if gca ~= udata.poly.THIS
  return;
end

coords(1,1) = getfield(get(gca,'CurrentPoint'),{3});
coords(1,2) = getfield(get(gca,'CurrentPoint'),{1});

udata.poly.x=cat(1,udata.poly.x,coords(2));
udata.poly.y=cat(1,udata.poly.y,coords(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove old 'polygon' from axis

delete(udata.poly.h);


switch get(gcf,'SelectionType')
 case 'normal'
  
  % add point to the 'polygon' and draw it
  
  for i=1:length(udata.a_h)
    subplot(udata.a_h(i));
    hold on;
    udata.poly.h(i) = plot(udata.poly.x,udata.poly.y,'black',...
			   'EraseMode','xor',...
			   'ButtonDownFcn',...
			   'vis_trajgui([],''click'')',...
			   'LineWidth',2);
  end
 case 'alt'
  
  % The polygon is ready.
  
  
  udata.poly.x=cat(1,udata.poly.x,udata.poly.x(1));
  udata.poly.y=cat(1,udata.poly.y,udata.poly.y(1));
  
  for i=1:length(udata.a_h)
    subplot(udata.a_h(i));
    udata.poly.h(i) = plot(udata.poly.x,udata.poly.y,'black',...
			   'EraseMode','xor',...
			   'ButtonDownFcn',...
			   'vis_trajgui([],''click'')',...
			   'LineWidth',2);
  end
  
  tmp=sort(repmat((1:udata.msize(1))',udata.msize(2),1));
  tmp(:,2)=repmat((1:udata.msize(2))',udata.msize(1),1);
  tmp2=tmp;
  if strcmp(udata.lattice,'hexa');
    t=find(~rem(tmp(:,1),2));
    tmp(t,2)=tmp(t,2)+0.5;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find the nodes that lie inside polygon and change coordinates to
  % linear indices.
  
  in = find(inpolygon(tmp(:,2),tmp(:,1),udata.poly.x,udata.poly.y));
  in = sub2ind(udata.msize,tmp2(in,1),tmp2(in,2));
  
  colors=get(udata.tmp_patch(1),'FaceVertexCData');
  colors=reshape(colors,[size(colors,1) 3]);
  tmp=ones(length(in),1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  % set the color of the nodes just selected, black.
  
  colors(in,:)=tmp*[0 0 0];
  set(udata.tmp_patch,'FaceVertexCData',colors);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % find the points mapping to the nodes from data 
  
  inds = [];
  for i=1:length(in)
    inds=[inds;find(in(i) == udata.bmus)];
  end
  
  %%%%%%%%%%%%%%%%%%%  
  % plot marks to data
  
  set(udata.fig1,'UserData',udata);
  plot2data(inds);
  udata=get(udata.fig1,'UserData');
  udata.new_marks=union(udata.new_marks,inds);
  delete(udata.poly.h);
  udata.poly.h=[];
  udata.poly.x=[];
  udata.poly.y=[];
  udata.poly.THIS=[];
end

set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_data()

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

data.points=[];
data.nodes=[];
k=1;

for i=1:length(udata.d_mark(1).h)
  data.points(i).inds=get(udata.d_mark(1).h(i),'XData');
  data.points(i).color=get(udata.d_mark(1).h(i),'Color');
end

color=get(udata.tmp_patch(1),'FaceVertexCData');
color=reshape(color,[size(color,1) 3]);

for i=1:size(color,1)
  if all(~isnan(color(i,:)))
    tmp.ind=i;
    tmp.color=color(i,:);
    data.nodes(k)=tmp;
    k=k+1;
  end
end

answer=inputdlg('Enter the name of the output variable:','',1);

if isempty(answer) | isempty(answer{1})
  msgbox('Output is not set to workspace.');
  return;
else
  assignin('base',answer{1},data);
  disp(sprintf('Struct is set to the workspace as ''%s''.',answer{1}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function load_data()

answer = inputdlg('Enter the name of the struct to be loaded:','',1);

if isempty(answer) | isempty(answer{1})
  msgbox('Data is not loaded.');
  return;
end

data=evalin('base',answer{1});

if ~isstruct(data)  
  errordlg('Input variable must be a struct.');
  return;
end

tmp1 = fieldnames(data);
tmp2 = {'nodes','points'};

for i=1:length(tmp1)
  for j=1:length(tmp2);
    if ~any(strcmp(tmp2{j},tmp1))
      errordlg('Wrong type of struct.');
      return;
    end
  end
end

if ~isempty(data.points)
  tmp1=fieldnames(data.points(1));
end
if ~isempty(data.nodes)
  tmp2=fieldnames(data.nodes(1));
end

for i=1:length(tmp1)
  if ~any(strcmp(tmp1{i},{'inds','color'}))
    errordlg('Wrong type of struct.');
    return;
  end
end

for i=1:length(tmp2)
  if ~any(strcmp(tmp2{i},{'ind','color'}))
    errordlg('Wrong type of struct.');
    return;
  end
end

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

clear_markers;
remove_traj;

old = gcf;

for i=1:length(data.points)
  for j=1:length(udata.h);
    set(0,'CurrentFigure',udata.fig1);
    subplot(udata.h(j));
    ydata=getfield(get(udata.d(j),'YData'),{data.points(i).inds});
    udata.d_mark(j).h=[udata.d_mark(j).h;...
                       plot(data.points(i).inds,ydata,...
			    'Color',data.points(i).color,...
			    'LineStyle','none',...
			    'Marker','o',...
			    'ButtonDownFcn',...
			    'vis_trajgui([],''line_down'')')];  
    if all(data.points(i).color == [0 0 0])
      udata.new_marks=unique([udata.new_marks; (data.points(i).inds)']);
    else
      udata.all_marks=unique([udata.all_marks; (data.points(i).inds)']);
    end
  end
  if ~isempty(udata.fig2)
    set(0,'CurrentFigure',udata.fig2);
    for j=1:length(udata.h2)
      subplot(udata.h2(j));
      ydata=getfield(get(udata.d2(j),'YData'),{data.points(i).inds});
      udata.d_mark2(j).h=[udata.d_mark2(j).h;...
                          plot(data.points(i).inds,ydata,...
			       'Color',data.points(i).color,...
			       'LineStyle','none',...
			       'Marker','o',...
			       'ButtonDownFcn',...
			       'vis_trajgui([],''line_down'')')];
    end
  end
end     


set(0,'CurrentFigure',udata.sM_h);
color=get(udata.tmp_patch(1),'FaceVertexCData');
color=reshape(color,[size(color,1) 3]);
for i=1:length(data.nodes)
  color(data.nodes(i).ind,:)=data.nodes(i).color;
end
for i=1:length(udata.tmp_patch);
  set(udata.tmp_patch(i),'FaceVertexCData',color);
end

set(0,'CurrentFigure',old);
set(udata.fig1,'UserData',udata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fuzzy_traj(trajStruct,arg); 
%function fuzzy_traj(sM_h,sM,sD,interval,arg)


if isempty(arg)
  
  if strcmp(trajStruct.lattice,'hexa')
    udata.lattice='hexa';
    udata.form=vis_patch('hexa');
  else
    data.lattice='rect';
    udata.form=vis_patch('rect');
  end
  
  % interval=[1 size(trajStruct.primary_data,1)];
  
  l=size(udata.form,1);
  dim = size(trajStruct.primary_data,2);
  udata.a_h=[findobj(get(trajStruct.figure,'Children'),'Tag','Uplane');...
             findobj(get(trajStruct.figure,'Children'),'Tag','Cplane')];
  
  udata.sM_h=trajStruct.figure;
  udata.msize=trajStruct.msize;
  
  %%%%%%%%%%%%%%%
  %
  % constructing patch that is drawn above every plane in map
  %
  
  nx = repmat(udata.form(:,1),1,prod(udata.msize));
  ny = repmat(udata.form(:,2),1,prod(udata.msize));
  
  x_c=reshape(repmat(1:udata.msize(2),l*udata.msize(1),1),l,prod(udata.msize));
  y_c=repmat(repmat(1:udata.msize(1),l,1),1,udata.msize(2));
  
  if strcmp(udata.lattice,'hexa')
    t = find(~rem(y_c(1,:),2));
    x_c(:,t)=x_c(:,t)+.5;
  end
  
  x_c=x_c+nx;
  y_c=y_c+ny;
  
  udata.orig_c=ones(prod(udata.msize),1)*[NaN NaN NaN];
  colors=reshape(udata.orig_c,[1 size(udata.orig_c,1) 3]);
  set(0,'CurrentFigure',trajStruct.figure);
  
  %%%%%%%%%
  % drawing 
  
  for i=1:length(udata.a_h);
    subplot(udata.a_h(i));
    v=caxis;
    udata.patch_h(i) =patch(x_c,y_c,colors,'EdgeColor','none');
    caxis(v);
  end
  
  
  udata.orig_x=get(udata.patch_h(1),'XData');
  udata.orig_y=get(udata.patch_h(1),'YData');
  
  %  if interval(1) < 1 | interval(2) > size(trajStruct.primary_data,1)
  %    error('Invalid argument ''interval''.');
  %  end
  
  x=1:size(trajStruct.primary_data,1);
  udata.fig1=figure;
  set(udata.fig1,'KeyPressFcn',...
                 'vis_trajgui([],''move_fuzzy'')');
  for i=1:size(trajStruct.primary_data,2)
    subplot(size(trajStruct.primary_data,2),1,i);
    udata.h(i)=gca;
    set(udata.h(1),'XTick',[]);
    udata.d(i)=plot(x,trajStruct.primary_data(:,i));
    l_x=1;
    lims(1) = round(min(trajStruct.primary_data(:,i)));
    lims(2) = round(max(trajStruct.primary_data(:,i)));
    udata.l(i) = line([l_x l_x],lims,'Color','red','EraseMode','xor');
  end
  
  udata.l_x = l_x;  
  udata.interval=[1 size(trajStruct.bmus,2)];
  
  tmp.fig1=udata.fig1;
  %  [K,P] = estimate_kernels(sM,sD);
  %  udata.K=K;
  %  udata.P=P;
  %  udata.codebook=sM.codebook;
  %  udata.data=sD.data;
  udata.bmus=trajStruct.bmus;
  set(udata.fig1,'UserData',udata);
  set(trajStruct.figure,'UserData',tmp);  
  draw_fuzzy(l_x);
  return;
end

move_fuzzy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_fuzzy(x)

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

inds = find(udata.bmus(:,x));
[row col] = ind2sub(udata.msize,inds);
if strcmp(udata.lattice,'hexa')
  t=find(~mod(row,2));
  col(t)=col(t)+0.5;
end

color=udata.orig_c;
color=reshape(color,[size(color,1) 3]);
xdata=udata.orig_x;
ydata=udata.orig_y;
tmp= ones(size(xdata(:,1),1),1)*udata.bmus(inds,x)';
color(inds,:) = ones(length(inds),1)*[0 0 0];
xdata(:,inds) = udata.form(:,1)*ones(1,length(inds)).*tmp+ones(6,1)*col';
ydata(:,inds) = udata.form(:,2)*ones(1,length(inds)).*tmp+ones(6,1)*row';

set(udata.patch_h,'FaceVertexCData',color,'XData',xdata,'YData',ydata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function move_fuzzy

% moves pointer lines and draws fuzzy response.

udata=get(gcf,'UserData');
udata=get(udata.fig1,'UserData');

switch get(gcf,'CurrentCharacter');
 case {'<','>'}
  key = get(gcf,'CurrentCharacter');
  if key == '>'
    if udata.l_x + 1 > udata.interval(2) 
      return;
    end
    l_x = udata.l_x + 1;
  else
    if udata.l_x - 1 < udata.interval(1)
      return;
    end
    l_x = udata.l_x - 1;
  end
  draw_fuzzy(l_x);
  set(udata.l,'XData',[l_x l_x]);
  udata.l_x=l_x;
  set(udata.fig1,'UserData',udata);
 otherwise
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



