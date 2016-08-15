type Params
    mu::Float64
    tout::Vector{Float64}
    yout::Vector{Vector{Float64}}
end

function newton(_n::Ptr{Cint}, _x::Ptr{Cdouble}, _y::Ptr{Cdouble}, _f::Ptr{Cdouble},
    _tnk::Ptr{Void})
    p = unsafe_pointer_to_objref(_tnk)
    mu = p.mu
    n = unsafe_load(_n, 1)
    y = unsafe_wrap(Array, _y, n, false)
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

function solout(nr::Ptr{Cint}, xold::Ptr{Cdouble}, x::Ptr{Cdouble},
    _y::Ptr{Cdouble}, n::Ptr{Cint}, con::Ptr{Cdouble}, icomp::Ptr{Cint},
    nd::Ptr{Cint}, _tnk::Ptr{Void}, irtrn::Ptr{Cint},
    xout::Ptr{Cdouble})
    p = unsafe_pointer_to_objref(_tnk)
    push!(p.tout, unsafe_load(x, 1))
    y = copy(unsafe_wrap(Array, _y, unsafe_load(n, 1), false))
    push!(p.yout, y)
    return nothing
end

cnewton = cfunction(newton, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Void}))
csolout = cfunction(solout, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
    Ptr{Cint}, Ptr{Void}, Ptr{Cint},
    Ptr{Cdouble}))

libdopri = Libdl.dlopen(Dopri.lib)
c_dopri5 = Libdl.dlsym(libdopri, :c_dopri5)
c_dop853 = Libdl.dlsym(libdopri, :c_dop853)
Libdl.dlclose(libdopri)

@testset "Low-Level" begin
    p = Params(398600.4415, Vector{Float64}(), Vector{Vector{Float64}}())
    s0 = [-1814.0, -3708.0, 5153.0, 6.512, -4.229, -0.744]
    y = copy(s0)
    xend = 5402.582703094263
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
    _n = Ref{Cint}(n)
    _x = Ref{Cdouble}(x)
    _xend = Ref{Cdouble}(xend)
    _itol = Ref{Cint}(itol)
    _iout = Ref{Cint}(iout)
    _lwork = Ref{Cint}(lwork)
    _liwork = Ref{Cint}(liwork)
    _idid = Ref{Cint}(idid)
    _tnk = pointer_from_objref(p)

    ccall(c_dopri5, Void, (Ref{Cint}, Ptr{Void}, Ref{Cdouble}, Ptr{Cdouble},
        Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint},
        Ptr{Void}, Ref{Cint}, Ptr{Cdouble}, Ref{Cint}, Ptr{Cint},
        Ref{Cint}, Ptr{Void}, Ref{Cint}),
        _n, cnewton, _x, y, _xend, rtol, atol,
        _itol, csolout, _iout, work, _lwork, iwork,
        _liwork, _tnk, _idid)

    # Copy results
    tj5 = copy(p.tout)
    yj5 = copy(p.yout)

    # Reinitialization
    empty!(p.tout)
    empty!(p.yout)
    y = copy(s0)
    x = 0.0
    work = zeros(Float64, lwork)
    iwork = zeros(Int32, liwork)
    _x = Ref{Cdouble}(x)
    _idid = Ref{Cint}(idid)

    ccall(c_dop853, Void, (Ref{Cint}, Ptr{Void}, Ref{Cdouble}, Ptr{Cdouble},
        Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint},
        Ptr{Void}, Ref{Cint}, Ptr{Cdouble}, Ref{Cint}, Ptr{Cint},
        Ref{Cint}, Ptr{Void}, Ref{Cint}),
        _n, cnewton, _x, y, _xend, rtol, atol,
        _itol, csolout, _iout, work, _lwork, iwork,
        _liwork, _tnk, _idid)

    tj8 = copy(p.tout)
    yj8 = copy(p.yout)

    @test length(tf5) == length(tj5)
    @test tf5 ≈ tj5
    for (a,b) in zip(yf5, yj5)
        @test a ≈ b
    end

    @test length(tf8) == length(tj8)
    @test tf8 ≈ tj8
    for (a,b) in zip(yf8, yj8)
        @test a ≈ b
    end
end
