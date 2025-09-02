include("./VDATqubithalf.jl")

data_dir="./gene_data/"
mkpath(data_dir)
plot_dir="./gene_plots"
mkpath(plot_dir)
# first step, compute fancyO ∼ ξ
for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        save_ξ_O(r,norb,data_dir;show_progress=false,Nθ=1000)
    end
end

# save_tildeO_with_derivative_and_Otp(norb,r; Δs=0.0001:0.0001:0.25)

# we have several intermediate step
# this one only computes the tildeO
for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        save_tildeO(norb,r,data_dir; Δs=0.001:0.0001:0.25)
    end
end


# now, generate with tildeO and derivatives. This will also generate Δc

for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        save_tildeO_with_derivative(norb,r,data_dir; Δs=0.001:0.0001:0.25)
    end
end


# we now also generate the data for two peak density of states. This will also generate Δc

for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        save_tildeO_with_derivative_and_Otp(norb,r,data_dir; Δs=0.001:0.0001:0.25)        
    end
end

# more precise way to compute Δc
for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        save_Δc_precise(norb,r,data_dir)
    end
end


for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        tag="bethe"
        saveData(solve_r_norb(r,norb,cal_Δ_bethe,cal_A_bethe,cal_K_bethe,data_dir;ainit=1.9,binit=0.1),filename_results(r,norb,tag,data_dir))
    end
end


for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        tag="bethe"
        print("norb $(norb) r $(r)\n")
        saveData(find_Uc_Δc_xx_ξc(r,norb,tag,data_dir),filename_Uc_Δc_xx_ξc(r,norb,tag,data_dir))
    end
end




# finally, we illustrate how to solve for a general density of state, using a set of discretized energy points
# we can use the previous code to generate an optimal discretization.

inf_fn=gene_spline_band("./es_files/es_inf.dat")
ϵs=gene_ϵs(inf_fn,0.5)[2]
cal_Δ_bethe_numeric,cal_A_bethe_numeric,cal_K_bethe_numeric=gene_cal_Δ_A_K_fn(ϵs)

for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        tag="bethe_numeric"
        saveData(solve_r_norb(r,norb,cal_Δ_bethe_numeric,cal_A_bethe_numeric,cal_K_bethe_numeric,data_dir;ainit=1.9,binit=0.1),filename_results(r,norb,tag,data_dir))
    end
end


# now, we check the two peak density of states, and compare with using O_tp

# The total kinetic energy is K0=-norb
ϵs=[1.0]
cal_Δ_two_peak,cal_A_two_peak,cal_K_two_peak=gene_cal_Δ_A_K_fn(ϵs)

for norb in 2:5
    for r in [0,0.05,0.1,0.25]
        print("norb $(norb) r $(r)\n")
        tag="two_peak"
        saveData(solve_r_norb(r,norb,cal_Δ_two_peak,cal_A_two_peak,cal_K_two_peak,data_dir;ainit=1.9,binit=0.1),filename_results(r,norb,tag,data_dir))
    end
end





