# Determine pure Fortran results
ref_file = normpath(joinpath(dirname(@__FILE__), "..", "deps", "ref.txt"))
results = open(readlines, ref_file, "r")
ind1 = 0
ind2 = 0
ind3 = 0
ind4 = 0
ind5 = 0
for (i, r) in enumerate(results)
    if ismatch(r"ind1", r)
        ind1 = i
    elseif ismatch(r"ind2", r)
        ind2 = i
    elseif ismatch(r"ind3", r)
        ind3 = i
    elseif ismatch(r"ind4", r)
        ind4 = i
    elseif ismatch(r"ind5", r)
        ind5 = i
    end
end
tf5 = Float64[]
yf5 = Array{Vector{Float64}}(0)
for i = 1:ind1-1
    v = float(split(results[i]))
    push!(tf5, v[1])
    push!(yf5, v[2:7])
end
tf8 = Float64[]
yf8 = Array{Vector{Float64}}(0)
for i = ind1+1:ind2-1
    v = float(split(results[i]))
    push!(tf8, v[1])
    push!(yf8, v[2:7])
end
tf5spc = Float64[]
yf5spc = Array{Vector{Float64}}(0)
for i = ind2+1:ind3-1
    v = float(split(results[i]))
    push!(tf5spc, v[1])
    push!(yf5spc, v[2:7])
end
tf8spc = Float64[]
yf8spc = Array{Vector{Float64}}(0)
for i = ind3+1:ind4-1
    v = float(split(results[i]))
    push!(tf8spc, v[1])
    push!(yf8spc, v[2:7])
end
tf5all = Float64[]
yf5all = Array{Vector{Float64}}(0)
for i = ind4+1:ind5-1
    v = float(split(results[i]))
    push!(tf5all, v[1])
    push!(yf5all, v[2:7])
end
tf8all = Float64[]
yf8all = Array{Vector{Float64}}(0)
for i = ind5+1:length(results)
    v = float(split(results[i]))
    push!(tf8all, v[1])
    push!(yf8all, v[2:7])
end
