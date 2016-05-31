using Compat

function getcompiler()
    which = is_unix() ? `which` : `where`
    if success(`$which ifort`)
        return "ifort"
    elseif success(`$which gfortran`)
        return "gfortran"
    else
        error("No compatible Fortran compiler found.")
    end
end

function unixbuild(compiler, path, ext)
    if compiler == "ifort"
        run(`ifort -O3 -xHost -ipo -fpic -c dopri.f90 dop853.f dopri5.f`)
        run(`ifort -dynamiclib -O3 -xHost -ipo -fpic -o $path/libdopri.$ext
            dopri.o dop853.o dopri5.o`)
        run(`ifort -o testrunner testrunner.f90 -ldopri -L$path`)
    elseif compiler == "gfortran"
        run(`gfortran -O3 -fpic -c dopri.f90 dop853.f dopri5.f`)
        run(`gfortran -shared -O3 -fpic -o $path/libdopri.$ext
            dopri.o dop853.o dopri5.o`)
        run(`gfortran -o testrunner -O3 testrunner.f90 -ldopri -L$path`)
    end
end

function windowsbuild(compiler, path, ext)
    error("Not yet supported.")
end

compiler = getcompiler()
path = splitdir(@__FILE__)[1]
ext = is_windows() ? "dll" : is_apple() ? "dylib" : "so"

build() = is_unix() ? unixbuild(compiler, path, ext) : windowsbuild(compiler, path, ext)
build()

for f in ["dopri.o", "dop853.o", "dopri5.o", "dopri.mod"]
    rm(f)
end
