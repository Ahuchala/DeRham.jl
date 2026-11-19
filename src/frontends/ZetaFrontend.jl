abstract type ZetaFrontend end

function polynomial end
function cohomology_quotient end
function weight end
function check_input end
function reduction_prescale end
function frobenius_series end
function frobenius_term_pole_order end

function reduction_prescale(P::CohomologyQuotient, m::Integer)
    m >= 1 || throw(ArgumentError("the pole order must be positive"))
    scale = QQ(1)
    for j = 1:(m-1)
        scale *= inv(reduction_scalar(P, j))
    end
    return scale
end
