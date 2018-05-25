function plot_vect(adjacentAssemblies,vect)

%% Read and create map
nass=length(adjacentAssemblies);
if nass~=length(adjacentAssemblies)
    error('ERROR : Wrong number of elements')
end
R=1;
r=R*sqrt(3)/2;

x=NaN*ones(nass,1);
y=NaN*ones(nass,1);

xdir=r*2*cos(pi/3.*(0:5));
ydir=r*2*sin(pi/3.*(0:5));

x(1)=0;
y(1)=0;

while sum(isnan(x))~=0
    Adj=adjacentAssemblies;
    for i=1:size(adjacentAssemblies,1)
        if ~isnan(x(i))
            for j=1:size(adjacentAssemblies,2)
                adj=Adj(i,j);
                if adj~=0
                    x(adj)=x(i)+xdir(j);
                    y(adj)=y(i)+ydir(j);
                    Adj(Adj==adj)=0;
                end
            end
        end
    end
end

%% plot results
%nw w sw se e ne
xhex=R*cos(pi/6+pi/3.*(0:5)); % x-coordinates of the vertices
yhex=R*sin(pi/6+pi/3.*(0:5)); % y-coordinates of the vertices

% Plotting part
axis equal
colormap(jet(1000))
c=ones(size(colormap));
c(0.25*length(c):0.65*length(c),:)=0;
colorbar
for i=1:length(x)
    col=round((length(c)-1)*(vect(i)-nanmin(vect))/(nanmax(vect)-nanmin(vect)));
    if col==0
        col=1;
    end
    if ~isnan(vect(i))
        p=patch(xhex+x(i),yhex+y(i),vect(i),'EdgeColor','none'); % make a hexagon at [2i,2j]
        tx=text(x(i), y(i), num2str(vect(i),3),'HorizontalAlignment','center','VerticalAlignment','middle','Color',c(col,:),'FontUnits', 'Normalized','FontSize',1.7e-2);
        hold on
    end
end
try
    filename=strfind(corename,'/');
    ax = gca;
    ax.Visible = 'off';
    title(corename(filename(end)+1:end),'Interpreter','none','Visible','on')
end