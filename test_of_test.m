it=[];
o=[];
F=cell(0);
nloc=100;
t=[];

while 1
    Fluxes_Input=randperm(nloc);
    tic
    test
    t(end+1,1)=toc;
    F{end+1,1}=Fluxes_Input;
    it(end+1,1)=output.iterations;
    o(end+1,1)=objval;
    save patterns F it o t
    pause(10)
end
    