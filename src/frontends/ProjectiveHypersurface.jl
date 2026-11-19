"""
    ProjectiveHypersurface(f)

The zeta frontend for the hypersurface defined by `f`.
"""
struct ProjectiveHypersurface{F} <: ZetaFrontend
    f::F
    P::HypersurfaceQuotient
end

function ProjectiveHypersurface(f)
    n = nvars(parent(f)) - 1
    d = total_degree(f)
    return ProjectiveHypersurface{typeof(f)}(f, HypersurfaceQuotient(n, d))
end

ProjectiveHypersurface{F}(f::F) where {F} = ProjectiveHypersurface(f)

polynomial(fe::ProjectiveHypersurface) = fe.f
cohomology_quotient(fe::ProjectiveHypersurface) = fe.P
weight(fe::ProjectiveHypersurface) = fe.P.n - 1

check_input(::ProjectiveHypersurface) = nothing

reduction_prescale(fe::ProjectiveHypersurface, m::Integer) =
    reduction_prescale(fe.P, m)

frobenius_series(::ProjectiveHypersurface, N, m) = computeD(N, m)

frobenius_term_pole_order(fe::ProjectiveHypersurface, m, j) =
    Int(characteristic(parent(polynomial(fe)))) * (m + j)
