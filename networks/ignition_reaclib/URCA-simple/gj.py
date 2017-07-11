from gauss_jordan import GaussJordan

GJ = GaussJordan(structure_file = 'sparsity.dat',
                 out_f95 = 'gauss_jordan_module.F90',
                 smp = False,
                 cse = False,
                 verbose = False)
