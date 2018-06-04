clear
addpath(genpath('/global/home/users/yvesrobert/CPLEX/cplex/matlab/x86-64_linux'));
addpath(genpath('Functions/'))
%% test
nloc=30;
nbatch=5;
ncycles=nbatch;
nass=nloc/ncycles;
Fluxes_Input=1:nloc;
Shuf_scheme=1:nbatch;
% Shuf_scheme=randperm(nbatch);
Batches=ceil([1:nloc]./nass);


%% Input
%Read Batch
validPositions = {[1], 1:6, 1:12, 1:18, 1:24, 1:30, 1:36, 1:42, [4:7, 11:15, 20:23, 27:31, 36:39, 43:47]}; %cell containing a vector for each ring. each vector contians the position numbers which are valid in each ring for assembly placement
pitch = 22.05; %cm
batchRadii = [0.0, (53.065+75.045)/2, (75.045+91.911)/2, (91.911+106.129)/2, (106.129+118.656)/2, (118.656+129.981)/2, (129.981+140.396)/2];
[~,Batches]=layoutOpt(validPositions,pitch,batchRadii);
%Read Power
assemblyPowerThreshold=1E4;
Shuf_scheme=[7 6 5 1 4 2 3];
nbatch=length(Shuf_scheme);
Q = readQ({'A_det0'});
lengthQ_original = length(Q);
[Q, map] = formMap(Q, assemblyPowerThreshold);
nloc=length(Q);
ncycles=nbatch;
nass=nloc/ncycles;
Fluxes_Input=Q./1e7;

%% Constraints
for k=1:ncycles
        Fluxes(:,k)=Fluxes_Input(Batches==Shuf_scheme(k));
end
npos=nass;

nvars=nass*ncycles + nass*npos*(ncycles-1) + 1;
%find how many constraints
nineqs = nass*(nass-1);
neqs = nass*ncycles+(ncycles-1)*nass + (ncycles-1)*npos;

ctype = '';
for i = 1:nvars
    if i <= nass*ncycles
        ctype(end+1) = 'C';
    elseif i <= nass*ncycles+nass*npos*(ncycles-1)
        ctype(end+1) = 'B';
    else
        ctype(end+1) = 'C';
    end
end

fprintf('\tBuilding constraint matrices\n');
idx_Aineq = 1; %keep track of the number of elements in Aineq matrix
idx_Aeq = 1;

%% Inequations
% Upper bound
constraint_idx = 1;
for i=1:nass %assembly
    for ip=1:nass
        if i~=ip
            for k=1:ncycles
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=(k-1)*nass+i;
                Aineq_v(idx_Aineq)=1;
                idx_Aineq = idx_Aineq+1;
                
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=(k-1)*nass+ip;
                Aineq_v(idx_Aineq)=-1;
                idx_Aineq = idx_Aineq+1;
            end
            Aineq_i(idx_Aineq)=constraint_idx;
            Aineq_j(idx_Aineq)=nass*ncycles + nass*npos*(ncycles-1) + 1;
            Aineq_v(idx_Aineq)=-1;
            idx_Aineq = idx_Aineq+1;
            
            bineq(constraint_idx) = 0;
            constraint_idx=constraint_idx+1;
        end
    end
end

%% Equalities
%fluence calculation
for i=1:nass
    Aeq_i(idx_Aeq) = i;
    Aeq_j(idx_Aeq) = i;
    Aeq_v(idx_Aeq) = 1;
    idx_Aeq = idx_Aeq+1;
    beq(i)=Fluxes(i,1);
end
for k=2:ncycles
    for i = 1:nass
        Aeq_i(idx_Aeq) = nass+(k-2)*nass+i;
        Aeq_j(idx_Aeq) = (k-1)*nass+i;
        Aeq_v(idx_Aeq) = 1;
        idx_Aeq = idx_Aeq+1;
        for j = 1:npos
            Aeq_i(idx_Aeq) = nass+(k-2)*nass+i;
            Aeq_j(idx_Aeq) = nass*ncycles+(k-2)*nass*npos+(i-1)*npos+j;
            Aeq_v(idx_Aeq) = -Fluxes(j,k);
            idx_Aeq = idx_Aeq+1;
        end
        beq(nass+(k-2)*nass+i) = 0;
    end
end

%exactly one location for an assembly per cycle
for k=2:ncycles
    for i = 1:nass
        for j = 1:npos
            Aeq_i(idx_Aeq) = ncycles*nass+(k-2)*nass+i;
            Aeq_j(idx_Aeq) = nass*ncycles+(k-2)*nass*npos+(i-1)*npos+j;
            Aeq_v(idx_Aeq) = 1;
            idx_Aeq = idx_Aeq+1;
        end
        beq(ncycles*nass+(k-2)*nass+i) = 1;
    end
end

%selection
for k=2:ncycles
    for j = 1:npos
        for i = 1:nass
            Aeq_i(idx_Aeq) = nass*ncycles+nass*(ncycles-1)+(k-2)*npos+j;
            Aeq_j(idx_Aeq) = nass*ncycles+(k-2)*nass*npos+(i-1)*npos+j;
            Aeq_v(idx_Aeq) = 1;
            idx_Aeq = idx_Aeq+1;
        end
        beq(nass*ncycles+nass*(ncycles-1)+(k-2)*npos+j) = 1;
    end
end

%% SOS
% sostype='';
% sosind=cell(1,nass*(ncycles-1));
% soswt=cell(1,nass*(ncycles-1));
% for k=2:ncycles
%     for i=1:nass
%         sostype(end+1)='1';
%         sosind{1,(k-2)*nass+i}=nass*ncycles + (k-2)*nass*npos + (i-1)*nass + [1:npos]';
%         soswt{1,(k-2)*nass+i}=Fluxes(:,k);
%     end
% end

%% Build and solve
fprintf('\tBuilding the constraint matrices\n');
%build sparse matrices for Aineq and Aeq
Aineq = sparse(Aineq_i, Aineq_j, Aineq_v, nineqs, nvars);
Aeq = sparse(Aeq_i, Aeq_j, Aeq_v, neqs, nvars);

%objective
c(nass*ncycles+nass*npos*(ncycles-1)+1)=1;

%print general problem parameters
fprintf('\t\tnumber of constraints = %i\n', neqs+nineqs);
fprintf('\t\t             equality = %i\n', neqs);
fprintf('\t\t           inequality = %i\n', nineqs);
fprintf('\t\tnumber of variables = %i\n', nvars);
fprintf('\t\t            integer = %i\n', sum(ctype == 'I'));
fprintf('\t\t             binary = %i\n', sum(ctype == 'B'));
fprintf('\t\t         continuous = %i\n', sum(ctype == 'C'));

%CPLEX solve
opts=cplexoptimset('display','on'); % Option to display iterations ('iter','on','off')
opts.exportmodel = 'model.lp'; % Name of the saved model
opts.workmem = 8000; % Maximum RAM to allocate for CPLEX
% [solutionvector, objval, status, output] = cplexmilp(c', Aineq, bineq', Aeq, beq', sostype, sosind, soswt, [], [], ctype,[],opts);
[solutionvector, objval, status, output] = cplexmilp(c', Aineq, bineq', Aeq, beq', [], [], [], [], [], ctype,[],opts);

fprintf('exit status = % i\n', status);
fprintf('solution time = % f\n', output.time);
fprintf('Objective value = % i\n', objval);

%% Post processing
Fluence(:,1)=Fluxes(:,1);
Acum_fluence=Fluence;
pos_in_batch(:,1)=1:npos;
for i=1:nass
    for k=2:ncycles
        Fluence(i,k)=solutionvector((k-1)*nass+i);
        pos_in_batch(i,k)=find(solutionvector(nass*ncycles+(k-2)*nass*npos+(i-1)*npos+(1:npos))==1);
        Acum_fluence(i,k)=Acum_fluence(i,k-1)+Fluence(i,k);
    end

end
pos=npos.*Shuf_scheme+pos_in_batch;
Max=solutionvector(end);
