# stat-em

run, in order, in a terminal:
\ generate_ground_truth.jl using_two second_var non_dim metrics
\ generate_emulator.jl using_two second_var non_dim metrics
\ generate_figures.jl using_two second_var non_dim metrics

where using_two toggles whether or not there's a second variable (tbh it's in some places hardcoded for two)
\second_var  is the CMIP name of the second variable (options for variables are tas, pr, huss, hurs)
\non_dim toggles normalization
\ metrics toggles metric terms in the basis and training

make sure project is pointing to the right place...
\ for example, I run:
\ julia --project=./StatEm generate_emulator.jl true hurs false false

which trains an emulator with tas and hurs and without normalization or metric terms