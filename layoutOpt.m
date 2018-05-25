function [ringPositionRadius,assemblyBatches]=layoutOpt(validPositions,pitch,batchRadii)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
%pre-optimization manipulation
%%%%%

numBatches = length(batchRadii);
numRings = length(validPositions);

%count the total number of assemblies specified
numAss = 0;
for ring = 1:length(validPositions)
    numAss = numAss + length(validPositions{ring});
end

numAssPerBatch = numAss/numBatches;

%make sure number of assemblies per batch is integral
if abs(mod(numAssPerBatch,1)) > 0.0001
    disp('error: number of assemblies per batch is not integral');
else
    numAssPerBatch = round(numAssPerBatch);
end

%initialize matrix to store assembly position (ring, position) and distance from core center
ringPositionRadius = cell(numAss,3);

%loop over each assembly specified and assign values for position and distance
assIdx = 1;
for ring = 1:length(validPositions)
    for pos = validPositions{ring}
        lpos = mod((pos-1),(ring-1)) + 1; %linear position, position of assembly along face of hex
        
        %assign the assembly type
        if lpos <= ring/2
            t = lpos;
        elseif lpos > ring/2
            t = ring-lpos+1;
        else
            disp('mistake');
        end
        
        %calculate distance from center of core and insert into cell
        ringPositionRadius{assIdx,1} = ring;
        ringPositionRadius{assIdx,2} = pos;
        ringPositionRadius{assIdx,3} = sqrt( ((ring-1)*pitch - (t-1)*pitch/2)^2 + ((t-1)*(3/2*pitch/sqrt(3)))^2 );
        
        assIdx = assIdx + 1;
    end
end

%create matrix with distance from each assembly to each r-z batch radius
distances = zeros(numAss, numBatches);
for ass = 1:numAss
    distances(ass,:) = abs(batchRadii-ringPositionRadius{ass,3});
end

%%%%%
%cplex stuff
%%%%%

%equality constraints
Aeq = zeros(numBatches+numAss, numAss*numBatches+6);
beq = ones(numBatches+numAss, 1);
for batch = 1:numBatches
    Aeq(batch, ((batch-1)*numAss+1):(batch*numAss)) = 1;
    beq(batch) = numAssPerBatch;
    for ass = 1:numAss
        Aeq(numBatches+ass, (batch-1)*numAss+ass) = 1;
    end
end

%inequality constraints
Aineq = zeros(6*2*numBatches, numAss*numBatches+6);
bineq = numAssPerBatch/6*ones(6*2*numBatches, 1);
bineq(2:2:end) = -bineq(2:2:end);
for batch = 1:numBatches
    for ass = 1:numAss
        ring = ringPositionRadius{ass, 1};
        pos = ringPositionRadius{ass, 2};
        if ring == 1
        else
            if pos <= (ring-1)*1
                Aineq((batch-1)*12+1, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+2, (batch-1)*numAss+ass) = -1;
                Aineq((batch-1)*12+11, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+12, (batch-1)*numAss+ass) = -1;
            elseif pos <= (ring-1)*2
                Aineq((batch-1)*12+1, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+2, (batch-1)*numAss+ass) = -1;
                Aineq((batch-1)*12+7, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+8, (batch-1)*numAss+ass) = -1;
            elseif pos <= (ring-1)*3
                Aineq((batch-1)*12+3, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+4, (batch-1)*numAss+ass) = -1;
                Aineq((batch-1)*12+7, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+8, (batch-1)*numAss+ass) = -1;
            elseif pos <= (ring-1)*4
                Aineq((batch-1)*12+3, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+4, (batch-1)*numAss+ass) = -1;
                Aineq((batch-1)*12+9, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+10, (batch-1)*numAss+ass) = -1;
            elseif pos <= (ring-1)*5
                Aineq((batch-1)*12+5, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+6, (batch-1)*numAss+ass) = -1;
                Aineq((batch-1)*12+9, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+10, (batch-1)*numAss+ass) = -1;
            elseif pos <= (ring-1)*6
                Aineq((batch-1)*12+5, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+6, (batch-1)*numAss+ass) = -1;
                Aineq((batch-1)*12+11, (batch-1)*numAss+ass) = 1;
                Aineq((batch-1)*12+12, (batch-1)*numAss+ass) = -1;
            else
                disp('error: something is wrong with the symmetry partitioning');
            end
        end
    end
end
Aineq(1:12:end, end-5) = -1;
Aineq(2:12:end, end-5) = -1;
Aineq(3:12:end, end-4) = -1;
Aineq(4:12:end, end-4) = -1;
Aineq(5:12:end, end-3) = -1;
Aineq(6:12:end, end-3) = -1;
Aineq(7:12:end, end-2) = -1;
Aineq(8:12:end, end-2) = -1;
Aineq(9:12:end, end-1) = -1;
Aineq(10:12:end, end-1) = -1;
Aineq(11:12:end, end) = -1;
Aineq(12:12:end, end) = -1;

%objective
c = zeros(numAss*numBatches+6, 1);
for batch = 1:numBatches
    c((batch-1)*numAss+1:batch*numAss) = ((numBatches+1)-batch)*distances(:,batch);
end
c(end-5:end) = 1;

%variable bounds
lb = -Inf*ones(numAss*numBatches+6,1);
ub = Inf*ones(numAss*numBatches+6,1);
%ub(end) = Inf;

%variable types
ctype = '';
for i = 1:numAss*numBatches
    ctype(end+1) = 'B';
end
ctype(end+1:end+6) = 'C';

%%%%%
%solve
%%%%%

[x, objval, status, output] = cplexmilp(c, Aineq, bineq, Aeq, beq, [], [], [], lb, ub, ctype);

%%%%%
%post-process
%%%%%

assemblyBatches = zeros(numAss, 1);
for batch  = 1:numBatches
    assemblyBatches = assemblyBatches + batch*x((batch-1)*numAss+1:batch*numAss);
end

%write plot file
fp = fopen([pwd '/Hexes_to_Plot.txt'],'w');
fprintf(fp, 'Ring, Position, Color, Value\n');
fprintf(fp, '#    #    #    #\n');

for ass = 1:numAss
    fprintf(fp, '%i %i SpringGreen %f\n', ringPositionRadius{ass,1}, ringPositionRadius{ass,2}, assemblyBatches(ass));
end

fclose(fp);