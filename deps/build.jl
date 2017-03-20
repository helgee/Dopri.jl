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
        run(`ifort -$(is_apple() ? "dynamiclib" : "shared") -O3 -xHost -ipo -fpic -o $path/libdopri.$ext dopri.f90`)
        run(`ifort -o testrunner testrunner.f90 -ldopri -L$path`)
    elseif compiler == "gfortran"
        run(`gfortran -shared -O3 -fpic -o $path/libdopri.$ext dopri.f90`)
        run(`gfortran -o testrunner -O3 testrunner.f90 -ldopri -L$path`)
    end
end

function windowsbuild(compiler, path, ext)
    if compiler == "gfortran"
        run(`gfortran -shared -O3 -fpic -o $path/libdopri.$ext dopri.f90`)
        run(`gfortran -o testrunner -O3 testrunner.f90 -ldopri -L$path`)
    end
end

compiler = getcompiler()
path = splitdir(@__FILE__)[1]
ext = is_windows() ? "dll" : is_apple() ? "dylib" : "so"

build() = is_unix() ? unixbuild(compiler, path, ext) : windowsbuild(compiler, path, ext)
build()

for f in ["dopri.mod"]
    rm(f)
end
