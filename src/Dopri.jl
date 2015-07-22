module Dopri

using Compat

export dop853, dopri5, dopricode

ext = @compat Dict(:Windows => "dll", :Darwin => "dylib", :Linux => "so")
const path = normpath(joinpath(splitdir(@__FILE__)[1],"..","deps"))
const lib = "$path/libdopri.$(ext[OS_NAME])"

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
    _nd::Ptr{Cint}, _rpar::Ptr{Cdouble}, _ipar::Ptr{Cint}, _irtrn::Ptr{Cint},
    _xout::Ptr{Cdouble})
    tnk = intstoobj(_ipar)
    if tnk.points != :last || tnk.S! != dummy
        n = unsafe_load(_n, 1)
        t = unsafe_load(_x, 1)
        told = unsafe_load(_xold, 1)
        y = unsafe_load(_y, n)
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

function _contd(i::Int, t, tnk::Thunk, _con::Ptr{Cdouble},
    _icomp::Ptr{Cint}, _nd::Ptr{Cint})
    if ~in(i, tnk.dense)
        error("No dense output available for element '$i'.")
    end
    return tnk.contd(i, @compat Float64(t), _con, _icomp, _nd)
end

dummy(xold, x, y, xout, irtrn, contd, params) = return nothing

function _fcn(_n::Ptr{Cint}, _x::Ptr{Cdouble}, _y::Ptr{Cdouble}, _f::Ptr{Cdouble},
    _rpar::Ptr{Cdouble}, _ipar::Ptr{Cint})
    tnk = intstoobj(_ipar)
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
            verbose::Bool=false)
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
            if ~verbose
                iwork[3] = -1
            end
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
            ipar = objtoints(tnk)

            ccall(($(sym), lib), Void, ($(doparg...),),
                &n, cfcn, &x, y, &xend, rtol, atol,
                &itol, csolout, &iout, work, &lwork, iwork,
                &liwork, rpar, ipar, &idid)

            if iout == 0
                push!(tout, tspan...)
                push!(yout, y0, y)
            end

            stats = @compat Dict{String,Int}("nfcn"=>iwork[17],"nstep"=>iwork[18],
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
    if WORD_SIZE == 64
        p = convert(Uint64, ptr)
        p1 = convert(Cuint, p >> 32)
        p2 = convert(Cuint, p << 32 >> 32)
        return [p1, p2]
    else
        return [convert(Uint32, ptr)]
    end
end

function intstoobj(arr)
    if WORD_SIZE == 64
        p1 = convert(Uint64, unsafe_load(arr, 1))
        p2 = convert(Uint64, unsigned(unsafe_load(arr,2)))
        ptr = convert(Ptr{Void}, p1 << 32 | p2)
    else
        ptr = convert(Ptr{Void}, unsafe_load(arr,1))
    end
    return unsafe_pointer_to_objref(ptr)
end

end # module
