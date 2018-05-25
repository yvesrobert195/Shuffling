fid=fopen('model.lp');
text = textscan(fid,'%s','delimiter','\n');
l=text{1};
for line=1:length(l)
    l{line}=[l{line} ' '];
    for k=1:ncycles
        l{line}=replace(l{line},['x' num2str(nass*ncycles+nass*ncycles*nloc+k) ' '],['D(' num2str(k) ') ']);
        for i=1:nass
            l{line}=replace(l{line},['x' num2str((k-1)*nass+i) ' '],['F(' num2str(i) ',' num2str(k) ') ']);
            for j=1:nloc
                l{line}=replace(l{line},['x' num2str(nass*ncycles+(k-1)*nass*nloc+(i-1)*nloc+j) ' '],['d(' num2str(i) ',' num2str(j) ',' num2str(k) ') ']);
            end
        end
    end
    l{line}=replace(l{line},['D(' num2str(nbatch) ') '],'D_tot');
end


fileID = fopen('model.txt','w');
fprintf(fileID,['Number of batches : ' num2str(nbatch) '\n']);
fprintf(fileID,['Number of locations/assemblies : ' num2str(nloc) '\n']);
fprintf(fileID,['Number of assemblies per batch : ' num2str(nass) '\n']);

for line = 1:length(l)
    fprintf(fileID,[l{line} '\n']);
end