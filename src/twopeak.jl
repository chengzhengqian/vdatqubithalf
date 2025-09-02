# functions for two-peak density of states, which is a convenient way to get U∼Δ

function filename_tildeO_with_derivatives_and_Otp(r,n_orb,data_dir)
    "$(data_dir)/Delta_tildeO_lnAbstildeO_xi_fancyFtilde2_fancyO_dtildeOdDelta_dtildeOdDeltaspline_dlnAbstildeOdDelta_dlnAbstildeOdDeltaspline_Otp_dOtpdDelta_dOtpdDeltaspline_n_orb_$(n_orb)_J_U_$(r).dat"    
end


# add the O from the two peak case, for an convenient check 
function save_tildeO_with_derivative_and_Otp(n_orb,r,data_dir; Δs=0.001:0.001:0.25)
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
    # add two peak case
    data_Otp=cal_Otp.(data_Δ,data_tildeO)
    # analytic derivative
    data_dOtpdΔ=cal_dOtpdDelta.(data_Δ,data_tildeO,data_dtildeOdΔ)
    # spline check
    Otp_Δ_fn=Spline1D(data_Δ,data_Otp)
    data_dOtpdΔspline=derivative(Otp_Δ_fn,data_Δ)
    # now, we save data in the follow 
    data_with_derivative_and_Otp=hcat(
        data_Δ,data_tildeO,data_lnAbstildeO,
        data_ξ,data_fancyFtilde2,data_fancyO,
        data_dtildeOdΔ,data_dtildeOdΔspline,
        data_dlnAbstildeOdΔ,data_dlnAbstildeOdΔspline,
        data_Otp,data_dOtpdΔ,data_dOtpdΔspline
    )
    saveData(data_with_derivative_and_Otp, filename_tildeO_with_derivatives_and_Otp(r,n_orb,data_dir))
    saveData(maximum([data_Δ[i] for i in 1:length(data_ξ) if (data_ξ[i]>1e-5)]), filename_Δc(r,n_orb,data_dir))
end




