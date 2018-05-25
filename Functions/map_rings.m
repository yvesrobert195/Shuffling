function ind=map_rings(adjacentAssemblies)
%% PART 1 : map concentric
R=1;
r=R*sqrt(3)/2;
xxdir=r*2*cos(pi/3.*(0:5));
yydir=r*2*sin(pi/3.*(0:5));
xxdir=[xxdir(3:end) xxdir(1:2)];
yydir=[yydir(3:end) yydir(1:2)];
%nw w sw se e ne

xxhexx=R*cos(pi/6+pi/3.*(0:5)); % xx-coordinates of the vertices
yyhexx=R*sin(pi/6+pi/3.*(0:5)); % yy-coordinates of the vertices

xx(1)=0;
yy(1)=0;
for i=2:20
    xx(end+1)=xx(1)+(i-1)*r*2;
    yy(end+1)=yy(1);
    for j=2:(i-1)*6
        xx(end+1)=xx(end)+xxdir(ceil((j-1)/(i-1)));
        yy(end+1)=yy(end)+yydir(ceil((j-1)/(i-1)));
    end
    i = i + 1;
end

%% PART 2 : map line by line
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
center=[(max(x)+min(x))/2,(max(y)+min(y))/2];
x=round(x-center(1),4);
y=round(y-center(2),4);

%% PART 3 : translation con <-> line
xx=round(xx',4);
yy=round(yy',4);
x=round(x,4);
y=round(y,4);
% 
for i=1:length(x)
    ind(i,1)=i;
    ind(i,2)=find(x(i)==xx & y(i)==yy);
end
