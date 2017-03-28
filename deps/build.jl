using Compat

function getcompiler()
    haskey(ENV, "FC") && return ENV["FC"]

    if success(`which ifort`)
        return "ifort"
    elseif success(`which gfortran`)
        return "gfortran"
    else
        error("No compatible Fortran compiler found.")
    end
end

function unixbuild(compiler, path, ext)
    if compiler == "ifort"
        run(`ifort -$(is_apple() ? "dynamiclib" : "shared") -O3 -xHost -ipo -fpic -o $path/libdopri.$ext dopri.f90`)
    elseif compiler == "gfortran"
        run(`gfortran -shared -O3 -fpic -o $path/libdopri.$ext dopri.f90`)
    end
end

function windowsbuild(path)
    # Compiled with GFortran from http://tdm-gcc.tdragon.net/
    # gfortran -static -O3 -fpic -m64 -o libdopri64.dll dopri.f90
    # gfortran -static -O3 -fpic -m32 -o libdopri32.dll dopri.f90
    download("https://dl.bintray.com/helgee/Dopri.jl/libdopri$(Sys.WORD_SIZE).dll",
        "$path/libdopri.dll")
end

compiler = is_unix() ? getcompiler() : ""
path = dirname(@__FILE__)
ext = is_apple() ? "dylib" : "so"

build() = is_unix() ? unixbuild(compiler, path, ext) : windowsbuild(path)
build()

for f in ["dopri.mod"]
    rm(f, force=true)
end
