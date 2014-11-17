module Dopri

export dop853, dopri5

ext = Dict(:Windows => "dll", :Darwin => "dylib", :Linux => "so")
const lib = "$(Pkg.dir())/Dopri/deps/libdopri.$(ext[OS_NAME])"

function solout(_nr::Ptr{Int32}, _xold::Ptr{Float64}, _x::Ptr{Float64},
    _y::Ptr{Float64}, _n::Ptr{Int32}, _con::Ptr{Float64}, _icomp::Ptr{Int32},
    _nd::Ptr{Int32}, _rpar::Ptr{Float64}, _ipar::Ptr{Int32}, _irtrn::Ptr{Int32},
    _xout::Ptr{Float64})
    tnk = intstoobj(unsafe_load(_ipar, 1), unsafe_load(_ipar, 2))
    n = unsafe_load(_n, 1)
    t = unsafe_load(_x, 1)
    push!(tnk.t, t)
    push!(tnk.y, copy(pointer_to_array(_y, n)))
    return nothing
end

function fcn(_n::Ptr{Int32}, _x::Ptr{Float64}, _y::Ptr{Float64}, _f::Ptr{Float64},
    _rpar::Ptr{Float64}, _ipar::Ptr{Int32})
    tnk = intstoobj(unsafe_load(_ipar, 1), unsafe_load(_ipar, 2))
    n = unsafe_load(_n, 1)
    t = unsafe_load(_x, 1)
    y = pointer_to_array(_y, n)
    f = pointer_to_array(_f, n)
    tnk.F!(f, t, y, tnk.params)
    return nothing
end

const fcnarg = (Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    Ptr{Float64}, Ptr{Int32})
const soloutarg = (Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
    Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32},
    Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
    Ptr{Float64})

cfcn = cfunction(fcn, Void, fcnarg)
csolout = cfunction(solout, Void, soloutarg)

type Thunk
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
    F!::Function
    params
end

syms = ["c_dopri5", "c_dop853"]
fcns = [:dopri5, :dop853]
for (fn, sym) in zip(fcns, syms)
    @eval begin
        function $(fn)(F!::Function, y0, tspan;
            params=[], atol=fill(sqrt(eps()), length(y0)),
            points=:all, rtol=fill(1e-6, length(y0)))
            x = tspan[1]
            xend = tspan[end]
            y = copy(y0)
            n = length(y0)
            tout = Float64[]
            yout = Array(Vector{Float64},0)
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
            if points == :all || points == :specified
                iout = 1
            elseif points == :last
                iout = 0
            else
                error("Unknown value for 'points'.")
            end
            lwork = 11*n + 8*n + 21
            liwork = n + 21
            work = zeros(Cdouble, lwork)
            iwork = zeros(Cint, liwork)
            rpar = zeros(Cdouble, 1)
            tnk = Thunk(tout, yout, F!, params)
            ipar = objtoints(tnk)
            idid = 0

            ccall(($(sym), lib), Void, (Ptr{Int32}, Ptr{Void}, Ptr{Float64},
                Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                Ptr{Int32}, Ptr{Void}, Ptr{Int32}, Ptr{Float64},
                Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
                Ptr{Cuint}, Ptr{Int32}), 
                &n, cfcn, &x, y, &xend, rtol, atol,
                &itol, csolout, &iout, work, &lwork, iwork,
                &liwork, rpar, ipar, &idid)

            if iout == 0
                push!(tout, x)
                push!(yout, y)
            end
            return tout, yout
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
