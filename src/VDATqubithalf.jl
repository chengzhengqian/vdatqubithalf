# solving the multi-orbital Hubbard model at half-filling using the Δ parametrization

using LinearAlgebra
using DelimitedFiles
using NLsolve
using Dierckx
using Optim
using QuadGK
using Arpack
using Statistics
using Roots    
# include necessary functions

# fancyFtilde^2 
include("./gene_code/fancyR.jl")
# The following three functions are for the half-filling cases.
include("./gene_code/Aflat.jl")
include("./gene_code/ximax.jl")
include("./gene_code/xigutz.jl")
include("./gene_code/dfancyRdDelta.jl")
# the density distribution
include("./gene_code/neabove.jl")
include("./gene_code/bethedos.jl")
# two peak case
include("./gene_code/Otp.jl")
include("./gene_code/dOtpdDelta.jl")
Power(x,y)=x^y
Sqrt=sqrt
Pi=pi
saveData(data,filename)= open(filename,"w") do io writedlm(io,data) end
loadData(filename)=readdlm(filename)
linspace(start,stop,length)=range(start,stop=stop,length=length)



include("./op.jl")

# mkdir(data_dir)
# we now set data_dir as explicit parameters.

function filename_ξ_O(r,n_orb,data_dir)
    "$(data_dir)/theta_xi_O_n_orb_$(n_orb)_J_U_$(r).dat"
end


function save_ξ_O(r,n_orb,data_dir;show_progress=false,Nθ=1000)
    data=cal_ξ_O(r,n_orb;show_progress=show_progress,Nθ=Nθ)
    saveData(data,filename_ξ_O(r,n_orb,data_dir))
end


include("./tildeO.jl")

function save_tildeO(n_orb,r,data_dir; Δs=0.001:0.001:0.25)
    fancyO_ξ_fn=get_fancyO(r,n_orb,data_dir)
    data=[[Δ,cal_tildeO_ξ_fancyFtilde2_fancyO(Δ,fancyO_ξ_fn)...] for Δ in Δs]
    saveData(data,"$(data_dir)/Delta_tildeO_xi_fancyFtilde2_fancyO_n_orb_$(n_orb)_J_U_$(r).dat")
end

include("./tildeOderivative.jl")

function filename_tildeO_with_derivatives(r,n_orb,data_dir)
    "$(data_dir)/Delta_tildeO_lnAbstildeO_xi_fancyFtilde2_fancyO_dtildeOdDelta_dtildeOdDeltaspline_dlnAbstildeOdDelta_dlnAbstildeOdDeltaspline_n_orb_$(n_orb)_J_U_$(r).dat"    
end

function filename_Δc(r,n_orb,data_dir)
    "$(data_dir)/Deltac_n_orb_$(n_orb)_J_U_$(r).dat"    
end


"""
save tildeO with derivative, computed in two different ways
also, for convenience, we also compute dln|tildeO|dΔ
we also save the deltac
"""
function save_tildeO_with_derivative(n_orb,r,data_dir; Δs=0.001:0.001:0.25)
    fancyO_ξ_fn=get_fancyO(r,n_orb,data_dir)
    data=[[Δ,cal_tildeO_ξ_fancyFtilde2_fancyO(Δ,fancyO_ξ_fn)...] for Δ in Δs]
    # analytic compute dtildeOdΔ
    data_Δ=[data_[1] for data_ in data]
    data_tildeO=[data_[2] for data_ in data]
    data_ξ=[data_[3] for data_ in data]
    data_fancyFtilde2=[data_[4] for data_ in data]
    data_fancyO=[data_[5] for data_ in data]
    data_dtildeOdΔ=cal_dtildeOdΔ.(data_Δ,data_ξ,data_fancyO)
    # use spline to compute dtildeOdΔ
    tildeO_Δ_fn=Spline1D(data_Δ,data_tildeO)
    data_dtildeOdΔspline=derivative(tildeO_Δ_fn,data_Δ)
    # we also compute the log version
    data_lnAbstildeO=log.(abs.(data_tildeO))
    data_dlnAbstildeOdΔ=data_dtildeOdΔ ./ data_tildeO
    # spline
    lnAbstildeO_Δ_fn=Spline1D(data_Δ,data_lnAbstildeO)
    data_dlnAbstildeOdΔspline=derivative(lnAbstildeO_Δ_fn,data_Δ)
    # now, we save data in the follow 
    data_with_derivative=hcat(
        data_Δ,data_tildeO,data_lnAbstildeO,
        data_ξ,data_fancyFtilde2,data_fancyO,
        data_dtildeOdΔ,data_dtildeOdΔspline,
        data_dlnAbstildeOdΔ,data_dlnAbstildeOdΔspline)
    saveData(data_with_derivative, filename_tildeO_with_derivatives(r,n_orb,data_dir))
    saveData(maximum([data_Δ[i] for i in 1:length(data_ξ) if (data_ξ[i]>1e-5)]), filename_Δc(r,n_orb,data_dir))
end


include("./twopeak.jl")


# we now try a more precise way to determine Δc (previously, just on the Δs grid)
function filename_Δc_precise(r,n_orb,data_dir)
    "$(data_dir)/Deltac_precise_n_orb_$(n_orb)_J_U_$(r).dat"    
end

"""
r=0.25
n_orb=2
"""
function save_Δc_precise(n_orb,r,data_dir)
    Δc_grid=loadData( filename_Δc(r,n_orb,data_dir))[1]
    data=loadData(filename_tildeO_with_derivatives(r,n_orb,data_dir))
    tildeO_metal=[data[i,2] for i in 1:size(data)[1] if data[i,1]<=Δc_grid][(end-1):end]
    Δ_metal=[data[i,1] for i in 1:size(data)[1] if data[i,1]<=Δc_grid][(end-1):end]
    tildeO_insulator=data[end,2]
    Δ_c=Δ_metal[2]+(tildeO_insulator-tildeO_metal[2])*(Δ_metal[2]-Δ_metal[1])/(tildeO_metal[2]-tildeO_metal[1])
    saveData(Δ_c,filename_Δc_precise(r,n_orb,data_dir))
end



# we now start to solve for a given density of states

include("bethe.jl")
include("general_dos.jl")

function solve_a_b(Δ,dlogtildeOdΔ,cal_Δ_fn,cal_A_fn;ainit=1.8,binit=0.1)
    function equations_ab!(F, x)
        a,b=x
        F[1] =(-8*a)-cal_A_fn(a,b)*b*dlogtildeOdΔ
        F[2] =cal_Δ_fn(a,b)-Δ
    end
    initial_guess = [ainit, binit]
    solution = nlsolve(equations_ab!, initial_guess)
    if(!converged(solution))
        print("not converged, with $(solution.zero)\n")
    end    
    solution.zero
end

# then we could compute other observables
function solve_a_b_A_U_O_E_K_Eloc(Δ,tildeO,dlogtildeOdΔ,cal_Δ_fn,cal_A_fn,cal_K_fn,norb;ainit=1.8,binit=0.1)
    a,b=solve_a_b(Δ,dlogtildeOdΔ,cal_Δ_fn,cal_A_fn;ainit=ainit,binit=binit)
    A=cal_A_fn(a,b)
    # Δ-cal_Δ_fn(a,b)
    K=cal_K_fn(a,b,norb)
    U=-norb*b/(2*A^3*tildeO)
    O=A^4*tildeO
    Eloc=U*O
    E=K+Eloc
    a,b,A,U,O,E,K,Eloc
end


function solve_r_norb(r,n_orb,cal_Δ_fn,cal_A_fn,cal_K_fn,data_dir;ainit=1.9,binit=0.1)
    # data=loadData("$(data_dir)/Delta_tildeO_xi_fancyR_barO_dtildeOdDeltaAnalytic_dtildeOdDeltaSpline_n_orb_$(n_orb)_J_U_$(r).dat")
    data=loadData(filename_tildeO_with_derivatives_and_Otp(r,n_orb,data_dir))
    data_Δ=data[:,1]
    data_tildeO=data[:,2]
    data_dlogtildeOdΔ=data[:,9]
    # data_dtildeOdΔ=data[:,7]
    # data_dlogtildeOdΔ .- data_dtildeOdΔ./data_tildeO
    a=ainit
    b=binit
    result=[]
    # idx=10
    for idx in 1:length(data_Δ)
        Δ=data_Δ[idx]
        dlogtildeOdΔ=data_dlogtildeOdΔ[idx]
        tildeO=data_tildeO[idx]
        a,b,A,U,O,E,K,Eloc=solve_a_b_A_U_O_E_K_Eloc(Δ,tildeO,dlogtildeOdΔ,cal_Δ_fn,cal_A_fn,cal_K_fn,n_orb;ainit=a,binit=b)
        print("Δ is $Δ and U is $U\n")
        push!(result,[Δ,a,b,A,U,O,E,K,Eloc])
    end
    result
end

"""
tag is the lattice name.
"""
function filename_results(r,n_orb,tag,data_dir)
    "$(data_dir)/Delta_a_b_A_U_O_E_K_Eloc_lattice_$(tag)_n_orb_$(n_orb)_J_U_$(r).dat"
end




# now, we retrieve various Δc and ξc (at Δc, independent of density of states)
# r=0.25
# norb=2
# tag="bethe"
# r=0.0
function find_Uc_Δc_xx_ξc(r,norb,tag,data_dir)
    # tag="bethe"
    data=loadData(filename_results(r,norb,tag,data_dir))
    # update to use more precise value of Δc, for the tildeO(Δ)
    Δc_MIT=loadData(filename_Δc(r,norb,data_dir))[1,1]
    Δc_MIT_precise=loadData(filename_Δc_precise(r,norb,data_dir))[1,1]
    data_metal_full=[(data[i,1],data[i,5]) for i in 1:size(data)[1] if data[i,1]<=Δc_MIT]
    Δc_Metal=argmax(x->x[2],data_metal_full)[1]
    if(Δc_Metal==Δc_MIT)
        Uc=data[findfirst(x->Δc_MIT==x,data[:,1]),5]
        #  in this case, ξc is zero
        ξc=0.0
        # Uc, Δc_MIT, Δc_Metal,Δc_Metal_1,Δc_Insulator_1
        # in this case, all critical Δc_XX are same
        return Uc,Δc_MIT_precise,Δc_MIT_precise,Δc_MIT_precise,Δc_MIT_precise,ξc
    end
    # avoid some initial small data
    data_metal_U=[data[i,5] for i in 1:size(data)[1] if data[i,1]<=Δc_Metal && data[i,1] >0.001]
    data_metal_E=[data[i,7] for i in 1:size(data)[1] if data[i,1]<=Δc_Metal  && data[i,1] >0.001]
    data_metal_Δ=[data[i,1] for i in 1:size(data)[1] if data[i,1]<=Δc_Metal  && data[i,1] >0.001]
    data_insulator_U=[data[i,5] for i in 1:size(data)[1] if data[i,1]>Δc_MIT]
    data_insulator_E=[data[i,7] for i in 1:size(data)[1] if data[i,1]>Δc_MIT]
    data_insulator_Δ=[data[i,1] for i in 1:size(data)[1] if data[i,1]>Δc_MIT]
    E_U_metal=Spline1D(data_metal_U,data_metal_E)
    E_U_insulator=Spline1D(data_insulator_U,data_insulator_E)
    Uc=find_zero(U->E_U_metal(U)-E_U_insulator(U),(data_insulator_U[1],data_metal_U[end]))
    Δ_U_metal=Spline1D(data_metal_U,data_metal_Δ)
    Δ_U_insulator=Spline1D(data_insulator_U,data_insulator_Δ)
    Δc_Metal_1=Δ_U_metal(Uc)
    Δc_Insulator_1=Δ_U_insulator(Uc)
    data_ξ=loadData(filename_tildeO_with_derivatives_and_Otp(r,norb,data_dir))[:,4]
    # data_Δ=loadData(filename_tildeO_with_derivatives_and_Otp(r,norb,data_dir))[:,1]
    idxc=maximum([i  for i in 1:length(data_ξ) if data_ξ[i]>1e-5])
    ξc=data_ξ[idxc]
    # Δc=data_Δ[idxc]
    Uc, Δc_MIT_precise, Δc_Metal,Δc_Metal_1,Δc_Insulator_1,ξc
end



# find_Uc_Δc_xx_ξc(0.0,norb,tag,data_dir)
# find_Uc_Δc_xx_ξc(0.25,norb,tag,data_dir)

function filename_Uc_Δc_xx_ξc(r,n_orb,tag,data_dir)
    "$(data_dir)/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_$(tag)_n_orb_$(n_orb)_J_U_$(r).dat"
end

