% Run and save data
FAV = linspace(0.1,0.5,41);  % invariant degree of polymerization
NQ=4;  % number of wavevector sets in calculating GAM4
NV=logspace(-1,4,15);

% write to file
filename='data/newgamdata';
if ~exist(filename,'file')
    outfile = fopen(filename, 'wt');
    for N=NV
        for FA=FAV
                 [GAM3,GAM4]=calcgamma(N,FA,NQ);
            fprintf(outfile,'%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
		                      N,FA,GAM3*N,GAM4(1:4)*N);
        end
    end
    fclose(outfile);
end