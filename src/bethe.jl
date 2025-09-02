# compute Δ,A, and K for the bethe lattice

function cal_Δ_bethe(a,b)
    integral, err = quadgk(x -> cal_bethedos(x)*cal_neabove(x,a,b), 0, 2, rtol=1e-10)
    integral
end

function cal_A_bethe(a,b)
    integral, err = quadgk(x -> cal_bethedos(x)*sqrt(cal_neabove(x,a,b)*(1.0-cal_neabove(x,a,b))), 0, 2, rtol=1e-10)
    2.0*integral
end

"""
compute the full kinetic energy (total kinetic energy per site)
cal_K_bethe(2.0,0.0,1)-(-8/3/pi)
"""
function cal_K_bethe(a,b,norb)
    integral, err = quadgk(x -> cal_bethedos(x)*x*(cal_neabove(x,a,b)-0.5), 0, 2, rtol=1e-10)
    integral*4*norb
end

