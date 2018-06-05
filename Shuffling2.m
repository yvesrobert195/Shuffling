clear
addpath(genpath('Functions/'))
%% test
nloc=12;
nbatch=3;
ncycles=nbatch;
nass=nloc/ncycles;
npos=nass;
Fluxes=1:nloc;
Shuf_scheme=randperm(nbatch);
Batches=ceil([1:nloc]./nass);


% Fluxes=randi(10,nloc,1);

%% Constraints
fprintf('\tInitialization of constraint matrices\n')

nvars=ncycles*nass+(ncycles-1)*npos*npos+1;
%find how many constraints
nineqs = nass*nass;
neqs = 0 ;

ctype = '';
for i = 1:nvars
    if i <= nass*ncycles
        ctype(end+1) = 'C';
    elseif i <= nass*ncycles+(ncycles-1)*npos*npos+1
        ctype(end+1) = 'B';
    else
        ctype(end+1) = 'C';
    end
end

fprintf('\tBuilding constraint matrices\n');
idx_Aineq = 1; %keep track of the number of elements in Aineq matrix
idx_Aeq = 1;

% Ineq
for i=1:nass %assembly
    for ip=1:nass
        for k=1:ncycles
            Aineq_i(idx_Aineq)=(i-1)*nass+ip;
            Aineq_j(idx_Aineq)=(k-1)*nass+i;
            if i~=ip
                Aineq_v(idx_Aineq)=1;
            else
                Aineq_v(idx_Aineq)=0;
            end
            idx_Aineq = idx_Aineq+1;
            
            Aineq_i(idx_Aineq)=(i-1)*nass+ip;
            Aineq_j(idx_Aineq)=(k-1)*nass+ip;
            if i~=ip
                Aineq_v(idx_Aineq)=-1;
            else
                Aineq_v(idx_Aineq)=0;
            end
            idx_Aineq = idx_Aineq+1;
        end
        Aineq_i(idx_Aineq)=(i-1)*nass+ip;
        Aineq_j(idx_Aineq)=ncycles*nass+(ncycles-1)*npos*npos+1;
        if i~=ip
            Aineq_v(idx_Aineq)=-1;
        else
            Aineq_v(idx_Aineq)=0;
        end
        idx_Aineq = idx_Aineq+1;
        bineq((i-1)*nass+ip) = 0;
    end
end

%Eq
for i=1:ncycles-1
    for j=1:nloc
        for k=1:loc
            Aeq_i(idx_Aeq)=(i-1)*nloc*nloc+(j-1)*nloc+k;
            Aeq_j(idx_Aeq)=(k-1)*nass+i;
            Aeq_v(idx_Aeq)=1;
        end
    end
end
          
%% Build and solve
fprintf('\tBuilding the constraint matrices\n');
%build sparse matrices for Aineq and Aeq
Aineq = sparse(Aineq_i, Aineq_j, Aineq_v, nineqs, nvars);

%objective
c(ncycles*nass+(ncycles-1)*npos*npos+1)=1;

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
[solutionvector, objval, status, output] = cplexmilp(c', Aineq, bineq', [], [], sostype, sosind, soswt, [], [], ctype,[],opts);

fprintf('exit status = % i\n', status);
fprintf('solution time = % f\n', output.time);
fprintf('Objective value = % i\n', objval);

read_model