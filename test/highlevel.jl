@testset "High-Level" begin
    # High-level interface
    function newton!(f, t, y, p)
        r = sqrt(y[1]*y[1]+y[2]*y[2]+y[3]*y[3])
        r3 = r*r*r
        f[1] = y[4]
        f[2] = y[5]
        f[3] = y[6]
        f[4] = -p.mu*y[1]/r3
        f[5] = -p.mu*y[2]/r3
        f[6] = -p.mu*y[3]/r3
    end

    type Params
        mu::Float64
        yout::Vector{Float64}
    end
    p = Params(mu, Float64[])
    tspan = [0.0, tp]
    tj5, yj5 = dopri5(newton!, s0, tspan, points=:last, params=p, maxstep=10.0)
    tj8, yj8 = dop853(newton!, s0, tspan, points=:last, params=p, maxstep=10.0)
    @test length(tj5) == 2
    @test length(yj5) == 2
    @test length(tj8) == 2
    @test length(tj8) == 2
    @test tj5 == tspan
    @test tj8 == tspan
    #= @test yj5[end] ≈ yf5[end] =#
    @test yj8[end] ≈ yf8[end]

    tj5, yj5 = dopri5(newton!, s0, tspan, points=:all, params=p)
    tj8, yj8 = dop853(newton!, s0, tspan, points=:all, params=p)
    @test length(tf5) == length(tj5)
    for (a,b) in zip(tf5, tj5)
        @test_approx_eq_eps a b tol
    end
    for (vf, vj) in zip(yf5, yj5)
        for (a, b) in zip(vf, vj)
            @test_approx_eq_eps a b tol
        end
    end

    @test length(tf8) == length(tj8)
    for (a,b) in zip(tf8, tj8)
        @test_approx_eq_eps a b tol
    end
    for (vf, vj) in zip(yf8, yj8)
        for (a, b) in zip(vf, vj)
            @test_approx_eq_eps a b tol
        end
    end

    # Test dense output
    tspan = collect(0.0:1.0:tp)
    tj5spc, yj5spc = dopri5(newton!, s0, tspan, points=:specified, params=p)
    tj8spc, yj8spc = dop853(newton!, s0, tspan, points=:specified, params=p)
    @test length(tf5spc) == length(tj5spc)
    for (a,b) in zip(tf5spc, tj5spc)
        @test_approx_eq_eps a b tol
    end
    for (vf, vj) in zip(yf5spc, yj5spc)
        for (a, b) in zip(vf, vj)
            @test_approx_eq_eps a b tol
        end
    end

    @test length(tf8spc) == length(tj8spc)
    for (a,b) in zip(tf8spc, tj8spc)
        @test_approx_eq_eps a b tol
    end
    for (vf, vj) in zip(yf8spc, yj8spc)
        for (a, b) in zip(vf, vj)
            @test_approx_eq_eps a b tol
        end
    end

    tspan = push!(collect(0.0:1.0:tp), tp)
    tj5all, yj5all = dopri5(newton!, s0, tspan, points=:all, params=p)
    tj8all, yj8all = dop853(newton!, s0, tspan, points=:all, params=p)
    @test length(tf5all) == length(tj5all)
    for (a,b) in zip(tf5all, tj5all)
        @test_approx_eq_eps a b tol
    end
    for (vf, vj) in zip(yf5all, yj5all)
        for (a, b) in zip(vf, vj)
            @test_approx_eq_eps a b tol
        end
    end

    @test length(tf8all) == length(tj8all)
    for (a,b) in zip(tf8all, tj8all)
        @test_approx_eq_eps a b tol
    end
    for (vf, vj) in zip(yf8all, yj8all)
        for (a, b) in zip(vf, vj)
            @test_approx_eq_eps a b tol
        end
    end

    # Test solout interface
    function solout!(told, t, y, contd, params)
        if told < 5000 < t
            push!(params.yout, 5000.0)
            for i = 1:6
                push!(params.yout, contd(i, 5000))
            end
            return dopricode[:abort]
        else
            return dopricode[:nominal]
        end
    end
    tj5, yj5 = dopri5(newton!, s0, tspan, solout=solout!, params=p)
    for (a,b) in zip([tf5spc[5001]; yf5spc[5001]], p.yout)
        @test_approx_eq_eps a b tol
    end
    p = Params(mu, Float64[])
    tj8, yj8 = dop853(newton!, s0, tspan, solout=solout!, params=p)
    for (a,b) in zip([tf8spc[5001]; yf8spc[5001]], p.yout)
        @test_approx_eq_eps a b tol
    end
end
