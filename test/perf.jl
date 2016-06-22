using Dopri
using BenchmarkTools

function newton!(f, t, y, mu)
    r = sqrt(y[1]*y[1]+y[2]*y[2]+y[3]*y[3])
    r3 = r*r*r
    f[1] = y[4]
    f[2] = y[5]
    f[3] = y[6]
    f[4] = -mu*y[1]/r3
    f[5] = -mu*y[2]/r3
    f[6] = -mu*y[3]/r3
end
tspan = [0.0, 5402.582703094263]
mu = 398600.4415
s0 = [-1814.0, -3708.0, 5153.0, 6.512, -4.229, -0.744]

b5 = @benchmark dopri5((f, t, y) -> newton!(f, t, y, mu), s0, tspan, points=:last)
b8 = @benchmark dop853((f, t, y) -> newton!(f, t, y, mu), s0, tspan, points=:last)
println("DOPRI5\n$b5")
println("DOP853\n$b8")

n = 1_00_000
@profile for i=1:n; dop853((f, t, y) -> newton!(f, t, y, mu), s0, tspan, points=:last); end
Profile.print()
