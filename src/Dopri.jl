module Dopri

export dop853, dopri5

ext = Dict(:Windows => "dll", :Darwin => "dylib", :Linux => "so")
const lib = "$(Pkg.dir())/Dopri/deps/libdopri.$(ext[OS_NAME])"

function _solout(_nr::Ptr{Cint}, _xold::Ptr{Cdouble}, _x::Ptr{Cdouble},
    _y::Ptr{Cdouble}, _n::Ptr{Cint}, _con::Ptr{Cdouble}, _icomp::Ptr{Cint},
    _nd::Ptr{Cint}, _rpar::Ptr{Cdouble}, _ipar::Ptr{Cint}, _irtrn::Ptr{Cint},
    _xout::Ptr{Cdouble})
    tnk = intstoobj(unsafe_load(_ipar, 1), unsafe_load(_ipar, 2))
    if tnk.points != :last
        n = unsafe_load(_n, 1)
        t = unsafe_load(_x, 1)
    end
    if tnk.points == :all
        push!(tnk.t, t)
        push!(tnk.y, copy(pointer_to_array(_y, n)))
    elseif tnk.points == :specified && t != 0.0
        told = unsafe_load(_xold, 1)
        times = tnk.t[(tnk.t .> told) & (tnk.t .<= t)]
        for ts in times
            if ts == t
                push!(tnk.y, copy(pointer_to_array(_y, n)))
            else
                yout = Float64[]
                for i=1:n
                    push!(yout, tnk.contd(i, ts, _con, _icomp, _nd))
                end
                push!(tnk.y, yout)
            end
        end
    end
    return nothing
end

function solout(xold, x, y, dense, params)
end

function _fcn(_n::Ptr{Cint}, _x::Ptr{Cdouble}, _y::Ptr{Cdouble}, _f::Ptr{Cdouble},
    _rpar::Ptr{Cdouble}, _ipar::Ptr{Cint})
    tnk = intstoobj(unsafe_load(_ipar, 1), unsafe_load(_ipar, 2))
    n = unsafe_load(_n, 1)
    t = unsafe_load(_x, 1)
    y = pointer_to_array(_y, n)
    f = pointer_to_array(_f, n)
    tnk.F!(f, t, y, tnk.params)
    return nothing
end

const fcnarg = (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cint})
const soloutarg = (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
    Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint},
    Ptr{Cdouble})
const doparg = (Ptr{Cint}, Ptr{Void}, Ptr{Cdouble},
    Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
    Ptr{Cint}, Ptr{Void}, Ptr{Cint}, Ptr{Cdouble},
    Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
    Ptr{Cuint}, Ptr{Cint})

cfcn = cfunction(_fcn, Void, fcnarg)
csolout = cfunction(_solout, Void, soloutarg)

type Thunk
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    F!::Function
    S!::Function
    points::Symbol
    params::Any
    contd::Function
end

syms = ["c_dopri5", "c_dop853"]
fcns = [:dopri5, :dop853]
dsyms = ["c_contd5", "c_contd8"]
dfcns = [:contd5, :contd8]
for (fn, sym, dfn, dsym) in zip(fcns, syms, dfcns, dsyms)
    @eval begin
        function $(fn)(F!::Function, y0, tspan;
            params::Any=[], atol::Vector{Float64}=fill(sqrt(eps()), length(y0)),
            points::Symbol=:all, rtol::Vector{Float64}=fill(1e-6, length(y0)),
            solout::Function=_solout)
            x = tspan[1]
            xend = tspan[end]
            y = copy(y0)
            n = length(y0)
            tout = Float64[]
            yout = Array(Vector{Float64},0)
            lwork = 11*n + 8*n + 21
            liwork = n + 21
            work = zeros(Cdouble, lwork)
            iwork = zeros(Cint, liwork)
            rpar = zeros(Cdouble, 1)
            idid = 0

            if length(rtol) != length(atol)
                error("Both 'atol' and 'rtol' must be provided in
                either scalar or vector form.")
            end
            if length(rtol) == 1
                rtol = [rtol]
                atol = [atol]
                itol = 0
            else
                itol = 1
            end
            if points == :all
                iout = 2
            elseif points == :specified
                iout = 2
                iwork[5] = n
                tout = tspan
                push!(yout, y0)
            elseif points == :last
                iout = 0
            else
                error("Unknown value for 'points'.")
            end
            tnk = Thunk(tout, yout, F!, solout, points, params, $(dfn))
            ipar = objtoints(tnk)

            ccall(($(sym), lib), Void, ($(doparg...),),
                &n, cfcn, &x, y, &xend, rtol, atol,
                &itol, csolout, &iout, work, &lwork, iwork,
                &liwork, rpar, ipar, &idid)

            if iout == 0
                push!(tout, tspan...)
                push!(yout, y0, y)
            end

            stats = Dict("nfcn"=>iwork[17],"nstep"=>iwork[18],
            "naccpt"=>iwork[19], "nrejct"=>iwork[20], "idid"=>idid)

            return tout, yout, stats
        end
        function $(dfn)(ii::Int, x::Float64, con::Ptr{Cdouble},
            icomp::Ptr{Cint}, nd::Ptr{Cint})
            ccall(($(dsym), lib), Cdouble,
                (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
                &ii, &x, con, icomp, nd)
        end
    end
end

function objtoints(obj)
    ptr = pointer_from_objref(obj)
    p = convert(Uint64, ptr)
    p1 = convert(Cuint, p >> 32)
    p2 = convert(Cuint, p << 32 >> 32)
    return [p1, p2]
end

function intstoobj(int1, int2)
    p1 = convert(Uint64, int1)
    p2 = convert(Uint64, unsigned(int2))
    ptr = convert(Ptr{Void}, p1 << 32 | p2)
    return unsafe_pointer_to_objref(ptr)
end

end # module
