__precompile__()

module Dopri

using Compat

export dop853, dopri5, dopricode

const lib = normpath(joinpath(dirname(@__FILE__), "..", "deps", "libdopri"))

function __init__()
    if Libdl.dlopen_e(lib) == C_NULL
        error("Please run Pkg.build(\"Dopri\").")
    end
    global const cfcn = cfunction(_fcn, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Void}))
    global const csolout = cfunction(_solout, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
        Ptr{Cint}, Ptr{Void}, Ptr{Cint},
        Ptr{Cdouble}))
end

immutable Irtrn
    value::Cint
end
dopricode = @compat Dict(:abort => Irtrn(-1),
    :altered => Irtrn(2),
    :nominal => Irtrn(0))

type DopriResults
    nfcn::Int
    nstep::Int
    naccpt::Int
    nrejct::Int
    idid::Int
end

type Thunk
    tspan::Vector{Float64}
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    F!::Function
    S!::Function
    points::Symbol
    contd::Function
    dense::Vector{Int}
    first::Bool
end

@compat abstract type DopriException <: Exception end

type DopriStiff <: DopriException
    msg::AbstractString
end

type DopriSmallStep <: DopriException
    msg::AbstractString
end

type DopriMaxStep <: DopriException
    msg::AbstractString
end

Base.show(io::IO, err::DopriException) = print(io, err.msg)

function _solout(_nr::Ptr{Cint}, _xold::Ptr{Cdouble}, _x::Ptr{Cdouble},
    _y::Ptr{Cdouble}, _n::Ptr{Cint}, _con::Ptr{Cdouble}, _icomp::Ptr{Cint},
    _nd::Ptr{Cint}, _tnk::Ptr{Void}, _irtrn::Ptr{Cint},
    _xout::Ptr{Cdouble})

    tnk = unsafe_pointer_to_objref(_tnk)::Thunk

    if tnk.points != :last || tnk.S! != dummy
        n = unsafe_load(_n, 1)
        t = unsafe_load(_x, 1)
        told = unsafe_load(_xold, 1)
        y = unsafe_wrap(Array, _y, n, false)
    end

    if tnk.points != :last && t != told
        if length(tnk.tspan) > 2
            times = tnk.tspan[(tnk.tspan .> told) .& (tnk.tspan .< t)]
            for ts in times
                push!(tnk.t, ts)
                yout = Array{Float64}(n)
                for i = 1:n
                    yout[i] = tnk.contd(i, ts, _con, _icomp, _nd)
                end
                push!(tnk.y, yout)
            end
        end
        if tnk.points == :all || t in tnk.tspan
            push!(tnk.t, t)
            push!(tnk.y, copy(y))
        end
    end
    if tnk.S! != dummy
        if t == told
            tnk.first = true
        else
            tnk.first = false
        end
        contd(i, t) = _contd(i, t, tnk, _con, _icomp, _nd)
        # Call the intermediate output function and assert that it
        # returns a valid return code.
        ret::Irtrn = tnk.S!(told, t, y, contd)
        unsafe_store!(_irtrn, ret.value, 1)
    end
    return nothing
end

function _contd(i::Int, _t, tnk::Thunk, _con::Ptr{Cdouble},
    _icomp::Ptr{Cint}, _nd::Ptr{Cint})
    if tnk.first
        error("Dense output function called at t=0.0.")
    end
    if ~in(i, tnk.dense)
        error("No dense output available for element '$i'.")
    end
    t = @compat Float64(_t)
    return tnk.contd(i, t, _con, _icomp, _nd)
end


dummy(xold, x, y, xout, irtrn, contd) = return nothing

function _fcn(_n::Ptr{Cint}, _x::Ptr{Cdouble}, _y::Ptr{Cdouble}, _f::Ptr{Cdouble},
    _tnk::Ptr{Void})
    tnk = unsafe_pointer_to_objref(_tnk)::Thunk
    n = unsafe_load(_n, 1)
    t = unsafe_load(_x, 1)
    y = unsafe_wrap(Array, _y, n, false)
    f = unsafe_wrap(Array, _f, n, false)
    tnk.F!(f, t, y)
    return nothing
end

syms = ["c_dopri5", "c_dop853"]
fcns = [:dopri5, :dop853]
dsyms = ["c_contd5", "c_contd8"]
dfcns = [:contd5, :contd8]
for (fn, sym, dfn, dsym) in zip(fcns, syms, dfcns, dsyms)
    @eval begin
        function $(fn)(F!::Function, y0, tspan;
            atol::Vector{Float64}=fill(sqrt(eps()), length(y0)),
            points::Symbol=:all, rtol::Vector{Float64}=fill(1e-6, length(y0)),
            solout::Function=dummy, dense::Vector{Int}=collect(1:length(y0)),
            verbose::Bool=false, safety::Float64=0.9, step_params::Vector{Float64}=[0.0,0.0],
            beta::Float64=0.0, maxstep::Real=tspan[end]-tspan[1], initstep::Float64=0.0,
            numstep::Int=100000, stiffness::Int=1000,
            )
            x = Float64(tspan[1])
            xend = Float64(tspan[end])
            y = Vector{Float64}(copy(y0))
            n = length(y)
            tout = Float64[x]
            yout = Array{Vector{Float64}}(0)
            push!(yout, copy(y))
            lwork = 11*n + 8*n + 21
            liwork = n + 21
            work = zeros(Float64, lwork)
            iwork = zeros(Int32, liwork)
            work[2] = safety
            work[3:4] = step_params
            work[5] = beta
            work[6] = maxstep
            work[7] = initstep
            if ~verbose
                iwork[3] = -1
            end
            iwork[1] = numstep
            iwork[4] = stiffness
            idid = 0

            if length(rtol) != length(atol)
                error("Both 'atol' and 'rtol' must be provided in
                either scalar or vector form.")
            end
            if length(rtol) == 1
                rtol = collect(rtol)
                atol = collect(atol)
                itol = 0
            else
                itol = 1
            end

            # iout=1 -> solout is called after every integration step
            # iout=2 -> solout is called and dense output is performed
            if points != :last || (solout != dummy && length(dense) != 0)
                iout = 2
                if points != :last
                    iwork[5] = n
                    dense = collect(1:n)
                elseif length(dense) != 0
                    iwork[5] = length(dense)
                    for (i, el) in enumerate(dense)
                        iwork[20+i] = el
                    end
                end
            elseif points == :last || (solout != dummy && length(dense) == 0)
                iout = 1
            elseif points == :last && solout == dummy
                iout = 0
            end
            tnk = Thunk(tspan, tout, yout, F!, solout, points, $(dfn), dense, false)

            _n = Ref{Cint}(n)
            _x = Ref{Cdouble}(x)
            _xend = Ref{Cdouble}(xend)
            _itol = Ref{Cint}(itol)
            _iout = Ref{Cint}(iout)
            _lwork = Ref{Cint}(lwork)
            _liwork = Ref{Cint}(liwork)
            _idid = Ref{Cint}(idid)
            _tnk = pointer_from_objref(tnk)
            ccall(($(sym), lib), Void, (Ref{Cint}, Ptr{Void}, Ref{Cdouble}, Ptr{Cdouble},
                Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint},
                Ptr{Void}, Ref{Cint}, Ptr{Cdouble}, Ref{Cint}, Ptr{Cint},
                Ref{Cint}, Ptr{Void}, Ref{Cint}),
                _n, cfcn, _x, y, _xend, rtol, atol,
                _itol, csolout, _iout, work, _lwork, iwork,
                _liwork, _tnk, _idid)

            if points == :last
                push!(tout, xend)
                push!(yout, y)
            end

            idid = _idid[]
            if idid == -1
                error("Input is not consistent.")
            elseif idid == -2
                throw(DopriMaxStep("Larger nmax needed."))
            elseif idid == -3
                throw(DopriSmallStep("Step size becomes too small."))
            elseif idid == -4
                throw(DopriStiff("Problem is probably stiff (interrupted)."))
            end

            stats = DopriResults(iwork[17], iwork[18], iwork[19], iwork[20], idid)

            return tout, yout, stats
        end

        function $(dfn)(ii::Int, x::Float64, con::Ptr{Cdouble},
            icomp::Ptr{Cint}, nd::Ptr{Cint})
            ccall(($(dsym), lib), Cdouble,
                (Ref{Cint}, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}),
                Ref{Cint}(ii), Ref{Cdouble}(x), con, icomp, nd)
        end
    end
end

end # module
