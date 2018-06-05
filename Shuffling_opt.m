clear
addpath(genpath('Functions/'))
%% test
nloc=90;
nbatch=10;
ncycles=nbatch;
nass=nloc/ncycles;
if mod(nass,1)~=0
    error('!!!')
end
Fluxes=1:nloc;
Shuf_scheme=randperm(nbatch);
Batches=ceil([1:nloc]./nass);


% Fluxes=randi(10,nloc,1);

%% Input
% %Read Batch
% validPositions = {[1], 1:6, 1:12, 1:18, 1:24, 1:30, 1:36, 1:42, [4:7, 11:15, 20:23, 27:31, 36:39, 43:47]}; %cell containing a vector for each ring. each vector contians the position numbers which are valid in each ring for assembly placement
% pitch = 22.05; %cm
% batchRadii = [0.0, (53.065+75.045)/2, (75.045+91.911)/2, (91.911+106.129)/2, (106.129+118.656)/2, (118.656+129.981)/2, (129.981+140.396)/2];
% [~,Batches]=layoutOpt(validPositions,pitch,batchRadii);
% %Read Power
% assemblyPowerThreshold=1E4;
% Shuf_scheme=[7 6 5 1 4 2 3];
% nbatch=length(Shuf_scheme);
% Q = readQ({'A_det0'});
% lengthQ_original = length(Q);
% [Q, map] = formMap(Q, assemblyPowerThreshold);
% nloc=length(Q);
% ncycles=nbatch;
% nass=nloc/ncycles;
% Fluxes=Q/sum(Q);
% %
% % Shuf_scheme=[12 11 10 9 8 1 7 6 2 3 4 5]';
% %

%% Constraints
fprintf('\tInitialization of constraint matrices\n')

nvars=nass*ncycles + nass*nloc*ncycles + ncycles-1 + 1;
%find how many constraints
nineqs = nass*(nass-1)*(ncycles-1)+nloc*ncycles + nass*(nass-1);
neqs = ncycles*nass + nass*ncycles;% + nass*ncycles ;

%initialize constraint matrices and vectors
% nelements_Aineq = 3*nass*nass*ncycles+nass*ncycles*nloc;
% Aineq_i = zeros(nelements_Aineq,1);
% Aineq_j = zeros(nelements_Aineq,1);
% Aineq_v = zeros(nelements_Aineq,1);
% bineq = zeros(nineqs, 1);
%
% nelements_Aeq = nass*nloc*sum(1:nbatch)+nass*ncycles*nloc + nass*ncycles*(nloc-nass);
% Aeq_i = zeros(nelements_Aeq,1);
% Aeq_j = zeros(nelements_Aeq,1);
% Aeq_v = zeros(nelements_Aeq,1);
% beq = zeros(neqs, 1);
% c = zeros(nvars,1);

ctype = '';
for i = 1:nvars
    if i <= nass*ncycles
        ctype(end+1) = 'C';
    elseif i <= nass*ncycles+nass*nloc*ncycles
        ctype(end+1) = 'B';
    elseif i <= nass*ncycles + nass*nloc*ncycles + ncycles-1
        ctype(end+1) = 'C';
    else
        ctype(end+1) = 'C';
    end
end

fprintf('\tBuilding constraint matrices\n');
idx_Aineq = 1; %keep track of the number of elements in Aineq matrix
idx_Aeq = 1;

% Ineq
constraint_idx = 1;
for k=1:ncycles-1
    for i=1:nass %assembly
        for ip=1:nass
            if i~=ip
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=k*nass+i;
                Aineq_v(idx_Aineq)=1;
                idx_Aineq = idx_Aineq+1;
                
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=(k-1)*nass+i;
                Aineq_v(idx_Aineq)=1;
                idx_Aineq = idx_Aineq+1;
                
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=k*nass+ip;
                Aineq_v(idx_Aineq)=-1;
                idx_Aineq = idx_Aineq+1;
                
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=(k-1)*nass+ip;
                Aineq_v(idx_Aineq)=-1;
                idx_Aineq = idx_Aineq+1;
                
                Aineq_i(idx_Aineq)=constraint_idx;
                Aineq_j(idx_Aineq)=nass*ncycles+nass*nloc*ncycles+k;
                Aineq_v(idx_Aineq)=-1;
                idx_Aineq = idx_Aineq+1;
                
                bineq(constraint_idx) = 0;
                constraint_idx=constraint_idx+1;
            end
        end
    end
end

%max one assembly per location
for k=1:ncycles
    for j = 1:nloc
        for i = 1:nass
            Aineq_i(idx_Aineq) = nass*(nass-1)*(ncycles-1)+(k-1)*nloc+j;
            Aineq_j(idx_Aineq) = nass*ncycles+(k-1)*nass*nloc+(i-1)*nloc+j;
            Aineq_v(idx_Aineq) = 1;
            idx_Aineq = idx_Aineq+1;
        end
        bineq(nass*(nass-1)*(ncycles-1)+(k-1)*nloc+j) = 1;
    end
end

% Upper bound
constraint_idx = nass*(nass-1)*(ncycles-1)+nloc*ncycles+1;
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
            Aineq_j(idx_Aineq)=nass*ncycles + nass*nloc*ncycles + ncycles-1 + 1;
            Aineq_v(idx_Aineq)=-1;
            idx_Aineq = idx_Aineq+1;
            
            bineq(constraint_idx) =0;
            constraint_idx=constraint_idx+1;
        end
    end
end

%% Eq
%fluence calculation
for k=1:ncycles
    for i = 1:nass
        Aeq_i(idx_Aeq) = (k-1)*nass+i;
        Aeq_j(idx_Aeq) = (k-1)*nass+i;
        Aeq_v(idx_Aeq) = 1;
        idx_Aeq = idx_Aeq+1;
        for j = 1:nloc
            Aeq_i(idx_Aeq) = (k-1)*nass+i;
            Aeq_j(idx_Aeq) = nass*ncycles+(k-1)*nass*nloc+(i-1)*nloc+j;
            if Batches(j)==Shuf_scheme(k)
                Aeq_v(idx_Aeq) = -Fluxes(j);
            else
                Aeq_v(idx_Aeq) =0;
            end
            idx_Aeq = idx_Aeq+1;
        end
        beq((k-1)*nass+i) = 0;
    end
end
%exactly one location for an assembly per cycle
for k=1:ncycles
    for i = 1:nass
        for j = 1:nloc
            Aeq_i(idx_Aeq) = ncycles*nass+(k-1)*nass+i;
            Aeq_j(idx_Aeq) = nass*ncycles+(k-1)*nass*nloc+(i-1)*nloc+j;
            Aeq_v(idx_Aeq) = 1;
            idx_Aeq = idx_Aeq+1;
        end
        beq(ncycles*nass+(k-1)*nass+i) = 1;
    end
end

% %selection
% for k=1:ncycles
%     for i=1:nass
%         for j=1:nloc
%             if Batches(j)~=Shuf_scheme(k)
%                 Aeq_i(idx_Aeq) = ncycles*nass+ncycles*nass+(k-1)*nass+i;
%                 Aeq_j(idx_Aeq) = nass*ncycles+(k-1)*nass*nloc+(i-1)*nloc+j;
%                 Aeq_v(idx_Aeq) = 1;
%                 idx_Aeq = idx_Aeq+1;
%             end
%         end
%         beq(ncycles*nass+ncycles*nass+(k-1)*nass+i) = 0;
%     end
% end

%% Build and solve
fprintf('\tBuilding the constraint matrices\n');
%build sparse matrices for Aineq and Aeq
Aineq = sparse(Aineq_i, Aineq_j, Aineq_v, nineqs, nvars);
Aeq = sparse(Aeq_i, Aeq_j, Aeq_v, neqs, nvars);

%objective
for k = 1:ncycles-1
    c(nass*ncycles+nass*nloc*ncycles+k) = 1;
end
c(nass*ncycles+nass*nloc*ncycles+ncycles-1+1)=1e6;

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
[solutionvector, objval, status, output] = cplexmilp(c', Aineq, bineq', Aeq, beq', [], [], [], [], [], ctype,[],opts);
% [solutionvector, objval, status, output] = cplexmilp(c', [], [], Aeq, beq', [], [], [], [], [], ctype,[],opts);

fprintf('exit status = % i\n', status);
fprintf('solution time = % f\n', output.time);
fprintf('Objective value = % i\n', objval);

%% Post-process
%read_model
if status==1
    for k=1:ncycles
        max_dif(k)=solutionvector(nass*ncycles+ncycles*nass*nloc+k);
        for i=1:nass
            for j=1:nloc
                deltas(i,j,k)=solutionvector(nass*ncycles+(k-1)*nass*nloc+(i-1)*nloc+j);
            end
            pos(i,k)=find(deltas(i,:,k)==1);
            batches(i,k)=Batches(pos(i,k));
        end
    end
    
    for i=1:nass
        for k=1:ncycles
            Fluence(i,k)=Fluxes(pos(i,k));
            Acum_Fluence(i,k)=solutionvector((k-1)*nass+i);
        end
        tot_Fluence(i,1)=sum(Fluence(i,:));
    end
else
    error('ERROR : NO SOLUTION FOUND')
end
% 
% read_model