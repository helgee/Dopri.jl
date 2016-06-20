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

function newton1!(f, t, y, mu)
    r = sqrt(y[1]*y[1]+y[2]*y[2]+y[3]*y[3])
    r3 = r*r*r
    f[1] = y[4]
    f[2] = y[5]
    f[3] = y[6]
    f[4] = -mu*y[1]/r3
    f[5] = -mu*y[2]/r3
    f[6] = -mu*y[3]/r3
end

type HighLevel
    mu::Float64
    yout::Vector{Float64}
end

function robertson!(f, t, y)
    f[1] = -0.04*y[1] + 10^4*y[2]*y[3]
    f[2] = 0.04*y[1] - 10^4*y[2]*y[3] - 3e7*y[2]^2
    f[3] = 3e7*y[2]^2
end

@testset "High-Level" begin
    mu = 398600.4415
    p = HighLevel(mu, Vector{Float64}())
    n1!(f, t, y) = newton1!(f, t, y, mu)
    s0 = [-1814.0, -3708.0, 5153.0, 6.512, -4.229, -0.744]
    tp = 5402.582703094263
    tspan = [0.0, tp]
    tj5, yj5 = dopri5(n1!, s0, tspan, points=:last)
    tj8, yj8 = dop853(n1!, s0, tspan, points=:last)
    @test length(tj5) == 2
    @test length(yj5) == 2
    @test length(tj8) == 2
    @test length(tj8) == 2
    @test tj5 == tspan
    @test tj8 == tspan
    @test yj5[end] ≈ yf5[end]
    @test yj8[end] ≈ yf8[end]

    n!(f, t, y) = newton!(f, t, y, p)
    tj5, yj5 = dopri5(n!, s0, tspan, points=:all)
    tj8, yj8 = dop853(n!, s0, tspan, points=:all)
    @test length(tf5) == length(tj5)
    @test tf5 ≈ tj5
    for (a, b) in zip(yf5, yj5)
        @test a ≈ b
    end

    @test length(tf8) == length(tj8)
    @test tf8 ≈ tj8
    for (a, b) in zip(yf8, yj8)
        @test a ≈ b
    end

    # Test dense output
    tspan = collect(0.0:1.0:tp)
    tj5spc, yj5spc = dopri5(n!, s0, tspan, points=:specified)
    tj8spc, yj8spc = dop853(n!, s0, tspan, points=:specified)
    @test length(tf5spc) == length(tj5spc)
    @test tf5spc ≈ tj5spc
    for (a, b) in zip(yf5spc, yj5spc)
        @test a ≈ b
    end

    @test length(tf8spc) == length(tj8spc)
    @test tf8spc ≈ tj8spc
    for (a, b) in zip(yf8spc, yj8spc)
        @test a ≈ b
    end

    tspan = push!(collect(0.0:1.0:tp), tp)
    tj5all, yj5all = dopri5(n!, s0, tspan, points=:all)
    tj8all, yj8all = dop853(n!, s0, tspan, points=:all)
    @test length(tf5all) == length(tj5all)
    @test tf5all ≈ tj5all
    for (a, b) in zip(yf5all, yj5all)
        @test a ≈ b
    end

    @test length(tf8all) == length(tj8all)
    @test tf8all ≈ tj8all
    for (a, b) in zip(yf8all, yj8all)
        @test a ≈ b
    end

    # Test solout interface
    function solout!(told, t, y, contd, params)
        if t == told
            return dopricode[:nominal]
        end
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
    s!(told, t, y, contd) = solout!(told, t, y, contd, p)
    tj5, yj5 = dopri5(n!, s0, tspan, solout=s!)
    @test [tf5spc[5001]; yf5spc[5001]] ≈ p.yout
    p = HighLevel(mu, Float64[])
    tj8, yj8 = dop853(n!, s0, tspan, solout=s!)
    @test [tf8spc[5001]; yf8spc[5001]] ≈ p.yout

    # Test error handling
    tspan = [0.0, tp]
    @test_throws ErrorException dopri5(n1!, s0, tspan, safety=-1.0)
    @test_throws Dopri.DopriMaxStep dopri5(n1!, s0, tspan, maxstep=1e-20)
    @test_throws Dopri.DopriMaxStep dopri5(n1!, s0, tspan, numstep=1)
    tspan = [0.0, 0.0]
    @test_throws Dopri.DopriSmallStep dopri5(n1!, s0, tspan)
    @test_throws Dopri.DopriStiff dopri5(robertson!, [1,0,0], [0,10^11])

    function solouterror!(told, t, y, contd)
        y = contd(1, t/2)
    end
    @test_throws ErrorException dopri5(n1!, s0, tspan, solout=solouterror!)
end
