# functions to compute dtildeOdΔ

# analytical formula to compute the derivative
function cal_dtildeOdΔ(Δ,ξ,fancyO)
    cal_dfancyRdDelta(Δ,ξ)*fancyO
end
