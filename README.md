# Qubit parametrization of VDAT for multi-orbital Hubbard model at half-filling

## Generate 1RDMF for a given J/U

### Include the file
https://github.com/chengzhengqian/vdatqubithalf/blob/36021b73a7c30b6ab8b3a08cadea901a9c2fbace/src/example.jl#L1

### Prepare the data directory 
https://github.com/chengzhengqian/vdatqubithalf/blob/36021b73a7c30b6ab8b3a08cadea901a9c2fbace/src/example.jl#L3-L4

### Generate $`\mathcal{O}(\xi)`$ for various $`N_{orb}`$ and  $`J/U`$
https://github.com/chengzhengqian/vdatqubithalf/blob/36021b73a7c30b6ab8b3a08cadea901a9c2fbace/src/example.jl#L8-L13

Results for $`\mathcal{O}(\xi)`$ 

![plot](./src/gene_plots/plot_xi_O_norb2_JU_0.0%200.05%200.25.png?raw=true)
![plot](./src/gene_plots/plot_xi_O_norb5_JU_0.0%200.05%200.25.png?raw=true)

### Compute $`\tilde{O}(\Delta)`$
https://github.com/chengzhengqian/vdatqubithalf/blob/36021b73a7c30b6ab8b3a08cadea901a9c2fbace/src/example.jl#L19-L24

Results
![plot](./src/gene_plots/plot_tildeO_Delta_norb2_JU_0.0%200.05%200.25.png?raw=true)

### Compute $`\frac{d\tilde{O}(\Delta)}{d\Delta}`$
https://github.com/chengzhengqian/vdatqubithalf/blob/36021b73a7c30b6ab8b3a08cadea901a9c2fbace/src/example.jl#L29-L34

Results for $`J/U=0.05`$ and $`N_{orb}=2`$
![plot](./src/gene_plots/plot_tildeO_and_derivative_in_log_norb2_JU_0.05.png?raw=true)

### Analyze the minimization of $`\xi`$ 
For $`J/U=0.05`$ and $`N_{orb}=2`$
![plot](./src/gene_plots/plot_fancyL_norb2_JU_0.05.png?raw=true)



## Solve for a given density of states

### Solve for Bethe lattice
https://github.com/chengzhengqian/vdatqubithalf/blob/36021b73a7c30b6ab8b3a08cadea901a9c2fbace/src/example.jl#L55-L61

Results
![plot](./src/gene_plots/plot_U_Delta_various_Norb_JU.png?raw=true)









