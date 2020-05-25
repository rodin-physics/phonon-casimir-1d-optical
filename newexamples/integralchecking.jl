using QuadGK
using Plots
pyplot()

const errtol = 1e-9;

# this is ω/Ω
k = 3
# a = 1

D1 = 3
D2 = 4

function A2(θ, k)
    return k^2 + sin(θ / 2)^2
end

function checkintegrand1(q2, θ, k, D2)
    Asq = A2(θ, k)
    I0 = exp(1im * q2 * D2)
    I3 = 1 / (Asq + sin(q2 / 2)  ^2)
    return I0 * I3
end

function checkIntegral1(θ, k, D2)
    f_check1(q2) = checkintegrand1(q2, θ, k, D2)
    res = quadgk(f_check1, 0, 2 * pi, atol = errtol)[1]
    return res
end

function asoln1(θ, k, D2)
    Asq = A2(θ, k)
    b = 2 * Asq + 1
    b2m1 = b ^ 2 - 1
    return 4 * pi * (b - sqrt(b2m1)) ^ D2 / sqrt(b2m1)
end


function checkintegrand2(q2, θ, k, D2)
    Asq = A2(θ, k)
    I0 = exp(1im * q2 * D2)
    I1 = sin(q2 / 2)  ^2
    I3 = 1 / (Asq + sin(q2 / 2)  ^2)
    return I0 * I1 * I3
end

function checkIntegral2(θ, k, D2)
    f_check2(q2) = checkintegrand2(q2, θ, k, D2)
    res = quadgk(f_check2, 0, 2 * pi, atol = errtol)[1]
    return res
end

function asoln2(θ, k, D2)
    Asq = A2(θ, k)
    b = 2 * Asq + 1
    b2m1 = b ^ 2 - 1
    rm = b - sqrt(b2m1)

    factor0 = pi
    factor1 = rm ^ D2
    factor2 = (2 - rm - 1 / rm)
    factor3 = 1 / sqrt(b2m1)
    return factor0 * factor1 * factor2 * factor3
end


# integrand1(-0.9878209999162516, k, D1, D2)
asoln1(-0.9878209999162516, k, D2)
checkIntegral1(-0.9878209999162516, k, D2)

r1 = map(q2 -> asoln1(q2, k, D2), 0:0.01: 2 * pi)
r2 = map(q2 -> checkIntegral1(q2, k, D2), 0:0.01: 2 * pi)

plot(0:0.01: 2 * pi, real.(r1))
plot(0:0.01: 2 * pi, real.(r2))


@time asoln2(-0.9878209999162516, k, D2)
@time checkIntegral2(-0.9878209999162516, k, D2)

r1p = map(q2 -> asoln2(q2, k, D2), 0:0.01: 2 * pi)
r2p = map(q2 -> checkIntegral2(q2, k, D2), 0:0.01: 2 * pi)

plot(0:0.01: 2 * pi, real.(r1p))
plot!(0:0.01: 2 * pi, real.(r2p))









Asq = A2(errtol, k)
I0 = 1 / (2 * pi)
I1 = exp(1im * errtol * D1) * (Asq - k)
I2 = (((Asq^2 + 1 / 2)^2 - 4) * (- 1im * 1im)) ^ (- 1/2)

function integral1(k, D1, D2)
    f_int1(θ) = integrand1(θ, k, D1, D2)
    res = quadgk(f_int1, 0, 2 * pi, atol = errtol)[1]
    return (res / (2 * pi))
end

function integral2(k, D1, D2)
    f_int2(θ) = integrand2(θ, k, D1, D2)
    res = quadgk(f_int2, - pi, pi, atol = errtol)[1]
    return (res / (2 * pi))
end

@time integral1(k, D1, D2)
@time integral2(k, D1, D2)
