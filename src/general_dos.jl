function gene_cal_Δ_A_K_fn(ϵs)
    function cal_Δ_fn(a,b)
        0.5*mean([cal_neabove(x,a,b)  for x in ϵs])
    end
    function cal_A_fn(a,b)
        mean([sqrt(cal_neabove(x,a,b)*(1.0-cal_neabove(x,a,b)))  for x in ϵs])
    end
    function cal_K_fn(a,b,norb)
        mean([x*(cal_neabove(x,a,b)-0.5)  for x in ϵs])*2*norb
    end
    cal_Δ_fn,cal_A_fn,cal_K_fn
end


# /home/chengzhengqian/share_workspace/czq_julia_package/VDATN3multi.jl/src/load_band.jl L1

function gene_spline_band(filename;scale=1.0)
    data=reshape(loadData(filename),:)*scale
    N_sample=size(data)[1]
    index=collect(linspace(0,1,N_sample))
    return Spline1D(index,data)
end

function gene_ϵs(e_fn,n;N_samples=40,N_minimal=4)
    index_points_below=max(trunc(Int64,n*N_samples),N_minimal)
    index_points_above=max(trunc(Int64,(1-n)*N_samples),N_minimal)
    index_below=collect(linspace(0,n,index_points_below+1))
    energy_below=[(integrate(e_fn,index_below[i],index_below[i+1]))/(index_below[i+1]-index_below[i]) for i in 1:index_points_below]
    index_above=collect(linspace(n,1,index_points_above+1))
    energy_above=[(integrate(e_fn,index_above[i],index_above[i+1]))/(index_above[i+1]-index_above[i]) for i in 1:index_points_above]
    [energy_below,energy_above]
end
