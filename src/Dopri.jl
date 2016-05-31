module Dopri

using Compat

export dop853, dopri5, dopricode

const path = normpath(joinpath(splitdir(@__FILE__)[1],"..","deps"))
const ext = is_windows() ? "dll" : is_apple() ? "dylib" : "so"
const lib = "$path/libdopri.$ext"

immutable Irtrn
    value::Cint
end
dopricode = @compat Dict(:abort => Irtrn(-1),
    :altered => Irtrn(2),
    :nominal => Irtrn(0))

type Thunk
    tspan::Vector{Float64}
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    F!::Function
    S!::Function
    points::Symbol
    params::Any
    contd::Function
    dense::Vector{Int}
end

function _solout(_nr::Ptr{Cint}, _xold::Ptr{Cdouble}, _x::Ptr{Cdouble},
    _y::Ptr{Cdouble}, _n::Ptr{Cint}, _con::Ptr{Cdouble}, _icomp::Ptr{Cint},
    _nd::Ptr{Cint}, _tnk::Ptr{Void}, _irtrn::Ptr{Cint},
    _xout::Ptr{Cdouble})

    tnk = unsafe_pointer_to_objref(_tnk)::Thunk

    if tnk.points != :last || tnk.S! != dummy
        n = unsafe_load(_n, 1)
        t = unsafe_load(_x, 1)
        told = unsafe_load(_xold, 1)
        y = pointer_to_array(_y, n)
    end

    if tnk.points != :last && t == told
        push!(tnk.t, t)
        push!(tnk.y, copy(pointer_to_array(_y, n)))
    elseif tnk.points != :last && t != told
        told = unsafe_load(_xold, 1)
        times = tnk.tspan[(tnk.tspan .> told) & (tnk.tspan .< t)]
        for ts in times
            push!(tnk.t, ts)
            yout = Float64[]
            for i=1:n
                push!(yout, tnk.contd(i, ts, _con, _icomp, _nd))
            end
            push!(tnk.y, yout)
        end
        if t in tnk.tspan || tnk.points == :all
            push!(tnk.t, t)
            push!(tnk.y, copy(pointer_to_array(_y, n)))
        end
    end
    if tnk.S! != dummy && t != told
        contd(i, t) = _contd(i, t, tnk, _con, _icomp, _nd)
        # Call the intermediate output function and assert that it
        # returns a valid return code.
        ret::Irtrn = tnk.S!(told, t, y, contd, tnk.params)
        unsafe_store!(_irtrn, ret.value, 1)
    end
    return nothing
end

function _contd(i::Int, _t, tnk::Thunk, _con::Ptr{Cdouble},
    _icomp::Ptr{Cint}, _nd::Ptr{Cint})
    if ~in(i, tnk.dense)
        error("No dense output available for element '$i'.")
    end
    t = @compat Float64(_t)
    return tnk.contd(i, t, _con, _icomp, _nd)
end

dummy(xold, x, y, xout, irtrn, contd, params) = return nothing

function _fcn(_n::Ptr{Cint}, _x::Ptr{Cdouble}, _y::Ptr{Cdouble}, _f::Ptr{Cdouble},
    _tnk::Ptr{Void})
    tnk = unsafe_pointer_to_objref(_tnk)::Thunk
    n = unsafe_load(_n, 1)
    t = unsafe_load(_x, 1)
    y = pointer_to_array(_y, n)
    f = pointer_to_array(_f, n)
    tnk.F!(f, t, y, tnk.params)
    return nothing
end

cfcn = cfunction(_fcn, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
Ptr{Void}))
csolout = cfunction(_solout, Void, (Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},
Ptr{Cint}, Ptr{Void}, Ptr{Cint},
Ptr{Cdouble}))

syms = ["c_dopri5", "c_dop853"]
fcns = [:dopri5, :dop853]
dsyms = ["c_contd5", "c_contd8"]
dfcns = [:contd5, :contd8]
for (fn, sym, dfn, dsym) in zip(fcns, syms, dfcns, dsyms)
    @eval begin
        function $(fn)(F!::Function, y0, tspan;
            params::Any=[], atol::Vector{Float64}=fill(sqrt(eps()), length(y0)),
            points::Symbol=:all, rtol::Vector{Float64}=fill(1e-6, length(y0)),
            solout::Function=dummy, dense::Vector{Int}=collect(1:length(y0)),
            verbose::Bool=false, safety::Float64=0.9, step_params::Vector{Float64}=[0.0,0.0],
            beta::Float64=0.0, maxstep::Float64=tspan[end]-tspan[1], initstep::Float64=0.0,
            numsteps::Int=100000, stiffness::Int=1000,
            )
            x = tspan[1]
            xend = tspan[end]
            y = copy(y0)
            n = length(y0)
            tout = Float64[]
            yout = Array(Vector{Float64},0)
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
            iwork[1] = numsteps
            iwork[4] = stiffness
            rpar = zeros(Cdouble, 1)
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
            elseif points == :last || solout == dummy
                iout = 0
            end
            tnk = Thunk(tspan, tout, yout, F!, solout, points, params, $(dfn), dense)

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
                push!(tout, tspan...)
                push!(yout, y0, y)
            end

            stats = @compat Dict{AbstractString,Int}("nfcn"=>iwork[17],"nstep"=>iwork[18],
            "naccpt"=>iwork[19], "nrejct"=>iwork[20], "idid"=>idid)

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
