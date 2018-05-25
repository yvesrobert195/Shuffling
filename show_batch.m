function lattice=show_batch(path)
file=fopen(path,'r');
lines= textscan(file,'%s','Delimiter','\n');
fclose(file);
lines=strtrim(lines{1});

i=1;
 while isempty(strfind(lines{i},'lat 1 '))
     i=i+1;
 end
 while isempty(strfind(lines{i},'0166'))
     i=i+1;
 end
 m=i;

while isempty(strfind(lines{m},'%%'))
    for j=1:8:length(lines{m})
            batch=lines{m}(j:j+3);
            lattice((m-i)/2+1,ceil(j/8))=str2double(batch(end-1:end));            
    end
    m=m+2;
end
while sum(~isnan(lattice(end,:)))==0
    lattice(end,:)=[];
end
while sum(~isnan(lattice(:,end)))==0
    lattice(:,end)=[];
end

% plot_core('',9,lattice)
% i=1;
%  while isempty(strfind(lines{i},'lat 1 '))
%      i=i+1;
%  end
%  while isempty(strfind(lines{i},'0166'))
%      i=i+1;
%  end
%  m=i;
% lattice=NaN*ones(ceil(length(lines{i})/8),3*ceil(length(lines{i})/8));
% while isempty(strfind(lines{m},'%%'))
%     for j=1:8:length(lines{m})
%             lattice((m-i)/2+1,(m-i)/2+2*ceil(j/8)-1)=str2double(lines{m}(j:j+3));            
%     end
%     m=m+2;
% end
% while sum(~isnan(lattice(end,:)))==0
%     lattice(end,:)=[];
% end
% while sum(~isnan(lattice(:,end)))==0
%     lattice(:,end)=[];
% end

% center=[ceil(size(lattice,1)/2),ceil(size(lattice,2)/2)];
% hexcor=[0 0;2 0;1 2; -1 2; -2 0; -1 -2;1 -2];

