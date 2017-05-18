addpath ../../../../matlab;
addpath ../../../../matlab/seq;
addpath ../../../../matlab/mike;
disp('Starting compilation...')
mcc -mC mutation_validator_wrapper
disp('Finished compilation.')
quit;
