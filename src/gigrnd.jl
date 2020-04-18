# Ported from PRScs/gigrnd.py

function psi(x, alpha, lam):
    f = -alpha*(cosh(x)-1)-lam*(exp(x)-x-1)
    return f
end

function dpsi(x, alpha, lam):
    f = -alpha*sinh(x)-lam*(exp(x)-1)
    return f


function g(x, sd, td, f1, f2):
    if (x >= -sd) && (x <= td)
        f = 1
    elseif x > td
        f = f1
    elseif x < -sd
        f = f2
    end

    return f
end


function gigrnd(p, a, b):
    # setup -- sample from the two-parameter version gig(lam,omega)
    p = float(p); a = float(a); b = float(b)
    lam = p
    omega = sqrt(a*b)

    if lam < 0
        lam = -lam
        swap = true
    else
        swap = false
    end

    alpha = sqrt(omega^2+lam^2)-lam

    # find t
    x = -psi(1, alpha, lam)
    if (x >= 1/2) && (x <= 2)
        t = 1
    elseif x > 2
        t = sqrt(2/(alpha+lam))
    elseif x < 1/2
        t = log(4/(alpha+2*lam))
    end

    # find s
    x = -psi(-1, alpha, lam)
    if (x >= 1/2) && (x <= 2)
        s = 1
    elseif x > 2
        s = sqrt(4/(alpha*cosh(1)+lam))
    elseif x < 1/2
        s = min(1/lam, log(1+1/alpha+sqrt(1/alpha^2+2/alpha)))
    end

    # find auxiliary parameters
    eta = -psi(t, alpha, lam)
    zeta = -dpsi(t, alpha, lam)
    theta = -psi(-s, alpha, lam)
    xi = dpsi(-s, alpha, lam)

    p = 1/xi
    r = 1/zeta

    td = t-r*eta
    sd = s-p*theta
    q = td+sd

    # random variate generation
    while true
        U = rand()
        V = rand()
        W = rand()
        if U < q/(p+q+r)
            rnd = -sd+q*V
        elseif U < (q+r)/(p+q+r)
            rnd = td-r*log(V)
        else
            rnd = -sd+p*log(V)
        end

        f1 = exp(-eta-zeta*(rnd-t))
        f2 = exp(-theta+xi*(rnd+s))
        if W*g(rnd, sd, td, f1, f2) <= exp(psi(rnd, alpha, lam))
            break
        end
    end

    # transform back to the three-parameter version gig(p,a,b)
    rnd = exp(rnd)*(lam/omega+sqrt(1+lam^2/omega^2))
    if swap
        rnd = 1/rnd
    end

    rnd = rnd/sqrt(a/b)
    return rnd
end
