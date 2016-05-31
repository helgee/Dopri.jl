@testset "Low-Level" begin
    # Low-level interface
    function newton(_n::Ptr{Int32}, _x::Ptr{Float64}, _y::Ptr{Float64}, _f::Ptr{Float64},
        _rpar::Ptr{Float64}, _ipar::Ptr{Int32})
        n = unsafe_load(_n, 1)
        y = pointer_to_array(_y, n)
        r = sqrt(y[1]*y[1]+y[2]*y[2]+y[3]*y[3])
        r3 = r*r*r
        unsafe_store!(_f, y[4], 1)
        unsafe_store!(_f, y[5], 2)
        unsafe_store!(_f, y[6], 3)
        unsafe_store!(_f, -mu*y[1]/r3, 4)
        unsafe_store!(_f, -mu*y[2]/r3, 5)
        unsafe_store!(_f, -mu*y[3]/r3, 6)
        return nothing
    end

    function solout(nr::Ptr{Int32}, xold::Ptr{Float64}, x::Ptr{Float64},
        _y::Ptr{Float64}, n::Ptr{Int32}, con::Ptr{Float64}, icomp::Ptr{Int32},
        nd::Ptr{Int32}, rpar::Ptr{Float64}, ipar::Ptr{Int32}, irtrn::Ptr{Int32},
        xout::Ptr{Float64})
        push!(tj, unsafe_load(x, 1))
        y = copy(pointer_to_array(_y, unsafe_load(n, 1)))
        push!(yj, y)
        return nothing
    end

    mu = 398600.4415
    s0 = [-1814.0, -3708.0, 5153.0, 6.512, -4.229, -0.744]
    y = copy(s0)
    tp = 5402.582703094263
    xend = tp
    x = 0.0
    idid = 0
    rtol = fill(1e-6, 6)
    atol = fill(1.4901161193847656e-8, 6)
    itol = 0
    iout = 1
    lwork = 200
    liwork = 100
    work = zeros(Float64, lwork)
    iwork = zeros(Int32, liwork)
    rpar = zeros(Float64, 1)
    ipar = zeros(Int32, 1)
    n = 6
    tj = Float64[]
    yj = Array(Vector{Float64},0)
    cnewton = cfunction(newton, Void, Dopri.fcnarg)
    csolout = cfunction(solout, Void, Dopri.soloutarg)

    ccall((:c_dopri5, Dopri.lib), Void, (Ptr{Int32}, Ptr{Void}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Void}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
        &n, cnewton, &x, y, &xend, rtol, atol,
        &itol, csolout, &iout, work, &lwork, iwork,
        &liwork, rpar, ipar, &idid)

    # Copy results
    tj5 = copy(tj)
    yj5 = copy(yj)

    # Reinitialization
    empty!(tj)
    empty!(yj)
    y = copy(s0)
    x = 0.0
    work = zeros(Float64, lwork)
    iwork = zeros(Int32, liwork)

    ccall((:c_dop853, Dopri.lib), Void, (Ptr{Int32}, Ptr{Void}, Ptr{Float64}, Ptr{Float64},
        Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
        Ptr{Void}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
        &n, cnewton, &x, y, &xend, rtol, atol,
        &itol, csolout, &iout, work, &lwork, iwork,
        &liwork, rpar, ipar, &idid)

    tj8 = copy(tj)
    yj8 = copy(yj)

    # The tolerance must not be too small due to the differences between
    # OPENLIBM, used by Julia,  and the LIBM (e.g. system or Intel) used by Fortran
    tol = 1e-5
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
end
