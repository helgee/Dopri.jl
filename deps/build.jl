function getcompiler()
    if success(`ifort --version`)
        return "ifort"
    elseif success(`gfortran --version`)
        return "gfortran"
    else
        error("No compatible Fortran compiler found.")
    end
end

path = "$(pwd())"
compiler = getcompiler()
ext = Dict(:Windows => "dll", :Darwin => "dylib", :Linux => "so")
dyn = Dict("ifort" => "-dynamiclib", "gfortran" => "-shared")
i8 = Dict("ifort" => "-i8", "gfortran" => "-fdefault-integer-8")
run(`$compiler -fpic -c dopri.f90 dop853.f dopri5.f`)
run(`$compiler $(dyn[compiler]) -fpic -o $path/libdopri.$(ext[OS_NAME])
    dopri.o dop853.o dopri5.o`)
run(`$compiler -o testrunner testrunner.f90 -ldopri -L$path`)

for f in ["dopri.o", "dop853.o", "dopri5.o", "dopri.mod"]
    rm(f)
end
