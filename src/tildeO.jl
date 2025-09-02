# functions to compute tildeO ∼ Δ


function get_fancyO(r,n_orb,data_dir)
    data=loadData(filename_ξ_O(r,n_orb,data_dir))
    ξ_data=data[end:-1:1,2]
    fancyO_data=data[end:-1:1,3]
    Spline1D(ξ_data,fancyO_data)
end


function cal_tildeO_ξ_fancyFtilde2_fancyO(Δ,fancyO_ξ_fn)
    ξmax=cal_ximax(Δ)
    result=optimize(ξ->cal_fancyR(Δ,ξ)*fancyO_ξ_fn(ξ),0.0,ξmax)
    tildeO=Optim.minimum(result)
    ξ=Optim.minimizer(result)
    # we should always compare with ξ=0
    # fancyR=fancyF^2 (removing the dependency of A)
    fancyR=cal_fancyR(Δ,ξ)
    fancyO=fancyO_ξ_fn(ξ)
    tildeOatomic=cal_fancyR(Δ,0)*fancyO_ξ_fn(0)
    if(tildeO<=tildeOatomic)
        return tildeO,ξ,fancyR,fancyO
    else
        return tildeOatomic,0.0,cal_fancyR(Δ,0),fancyO_ξ_fn(0)
    end    
end


