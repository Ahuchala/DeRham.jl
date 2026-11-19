@testset "depthfirst reduction policy" begin
    geometry = DeRham.ReductionEngineGeometry(5, 2)
    directions = DeRham.ReductionDirectionSet(
        (u, k, v) ->
            length(u) == length(v) && sum(v) == 3 && all(v .>= 0) && all(u .>= k .* v),
    )
    engine = DeRham.AccountingReductionEngine(geometry, directions)
    policy = DeRham.DepthFirstPolicy(3, Int[])
    inputs = [(key = ([4, 4, 4], 4), terms = []), (key = ([3, 4, 4], 4), terms = [])]

    terminal_chains =
        DeRham.run_reduction_policy!(policy, engine, [(4, inputs)], directions)

    result = DeRham.accounting_result(engine)
    @test length(terminal_chains) == 2
    @test all(chain -> chain.order == geometry.terminal_order, terminal_chains)
    @test result.live_chains == 0
    @test result.peak_live_chains == 1
    @test result.steps == 4
    @test result.chunks == 4
    @test all(
        event -> directions(event.key_before.u, event.chunk_length, event.direction),
        filter(event -> event.operation == :reduce, engine.events),
    )
end
