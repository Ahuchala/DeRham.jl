
function reducetransform(
    FT,
    N_m,
    S,
    f,
    P,
    fe::ZetaFrontend,
    pseudoInverseMat,
    p,
    params,
    cache,
    context,
)
    if params.algorithm == :pchunk
        reducetransform_pchunk(FT, N_m, S, f, P, pseudoInverseMat, p, cache, params, fe)
    elseif params.algorithm == :depthfirst
        reducetransform_depthfirst(
            FT,
            N_m,
            S,
            f,
            P,
            fe,
            pseudoInverseMat,
            p,
            cache,
            params,
            context,
        )
    elseif params.algorithm == :varbyvar
        reducetransform_varbyvar(
            FT,
            N_m,
            S,
            f,
            P,
            pseudoInverseMat,
            p,
            cache,
            params,
            context,
            fe,
        )
    elseif params.algorithm == :akr
        reducetransform_akr(FT, N_m, S, f, pseudoInverseMat, p, cache, params)
    else
        println(params.algorithm)
        throw(ArgumentError("Unsupported Algorithm: "))
    end
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    P::HypersurfaceQuotient,
    pseudoInverseMat,
    p,
    params,
    cache,
    context,
)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        ProjectiveHypersurface(f),
        pseudoInverseMat,
        p,
        params,
        cache,
        context,
    )
end

function reducetransform(
    FT,
    N_m,
    S,
    f,
    P::DoubleCoverQuotient,
    pseudoInverseMat,
    p,
    params,
    cache,
    context,
)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        DoubleCover(f),
        pseudoInverseMat,
        p,
        params,
        cache,
        context,
    )
end

function reducetransform(FT, N_m, S, f, pseudoInverseMat, p, params, cache, context)
    fe = ProjectiveHypersurface(f)
    P = cohomology_quotient(fe)
    return reducetransform(
        FT,
        N_m,
        S,
        f,
        P,
        fe,
        pseudoInverseMat,
        p,
        params,
        cache,
        context,
    )
end

function reducetransform(FT, N_m, S, f, pseudoInverseMat, p, params, cache)
    return reducetransform(FT, N_m, S, f, pseudoInverseMat, p, params, cache, nothing)
end
