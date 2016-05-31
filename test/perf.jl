using Dopri

function runner(f, n)
    tau = 0.0
    gc_enable(false)
    for i=1:n
        tau += @elapsed(f())
    end
    gc_enable(true)
    return tau/n
end

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

n = 10000
dopri5(newton!, s0, tspan, points=:all, params=mu)
dop853(newton!, s0, tspan, points=:all, params=mu)
tau = runner(()->dopri5(newton!, s0, tspan, points=:last, params=mu), n)
println("DOPRI5:\t$(tau) seconds per loop")
tau = runner(()->dop853(newton!, s0, tspan, points=:last, params=mu), n)
println("DOP853:\t$(tau) seconds per loop")
