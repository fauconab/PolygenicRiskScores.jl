# Ported from PRCcs/src/mcmc_gtb.py

backend(::Type{Array}) = KernelAbstractions.CPU()
mzeros(::Type{Array{T,N}}, args...) where {T,N} = zeros(T, args...)
mones(::Type{Array{T,N}}, args...) where {T,N} = ones(T, args...)
mrand(::Type{Array{T,N}}, args...) where {T,N} = rand(T, args...)
mrandn(::Type{Array{T,N}}, args...) where {T,N} = randn(T, args...)

function mcmc(a, b, phi, sst_df, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, beta_std, seed; verbose=false, Tval=Float64, Tarr=Array)
    Vec = Tarr{Tval, 1}
    Mat = Tarr{Tval, 2}

    # seed
    if seed !== nothing
        Random.seed!(seed)
    end

    # derived stats
    beta_mrg = convert(Vec, copy(sst_df.BETA))
    maf = convert(Vec, copy(sst_df.MAF))
    n_pst = (n_iter-n_burnin)/thin
    p = length(sst_df.SNP)
    n_blk = length(ld_blk)

    # initialization
    beta = mzeros(Vec, p)
    psi = mones(Vec, p)
    sigma = Tval(1.0)
    if phi === nothing
        phi = Tval(1.0)
        phi_updt = true
    else
        phi_updt = false
    end

    beta_est = mzeros(Vec, p)
    psi_est = mzeros(Vec, p)
    sigma_est = Tval(0.0)
    phi_est = Tval(0.0)

    a = Tval(a)
    b = Tval(b)

    Rlen = 2^2

    ld_blk = convert.(Ref(Mat), ld_blk)

    # MCMC
    for itr in 1:n_iter
        if itr % 100 == 0
            verbose && @info "(Chromosome $chrom) MCMC iteration $itr"
        end

        mm = 1; quad = 0.0
        for kk in 1:n_blk
            if blk_size[kk] == 0
                continue
            else
                idx_blk = mm:(mm+blk_size[kk]-1)
                dinvt = ld_blk[kk] .+ Diagonal(Tval(1.0) ./ psi[idx_blk])
                dinvt_chol = cholesky(dinvt).U
                beta_tmp = (transpose(dinvt_chol) \ beta_mrg[idx_blk]) .+ sqrt(sigma/n) .* mrandn(Vec, length(idx_blk))
                beta[idx_blk] = dinvt_chol \ beta_tmp
                quad += dot(transpose(beta[idx_blk]) * dinvt, beta[idx_blk])
                mm += blk_size[kk]
            end
        end

        err = max(n/2.0*(1.0-2.0*sum(beta.*beta_mrg)+quad), n/2.0*sum(beta .^ 2 ./ psi))
        sigma = 1.0/rand(Gamma((n+p)/2.0, 1.0/err))

        delta = convert(Vec, rand.(Gamma.(a+b, Array(1.0 ./ (psi .+ phi)))))

        R = mrand(Mat, length(psi), Rlen)
        gigrnd_F = gigrnd(backend(Tarr), length(psi))
        wait(gigrnd_F(psi, a-0.5, 2.0 .* delta, n .* (beta .^ 2) ./ sigma, R; ndrange=length(psi)))
        psi[psi .> 1] .= 1.0

        if phi_updt
            w = rand(Gamma(1.0, 1.0/(phi+1.0)))
            phi = rand(Gamma(p*b+0.5, 1.0/(sum(delta)+w)))
        end

        # posterior
        if (itr>n_burnin) && (itr % thin == 0)
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst
        end
    end

    # convert standardized beta to per-allele beta
    if !beta_std
        beta_est ./= sqrt.(2.0 .* maf .* (1.0 .- maf))
    end

    # print estimated phi
    if phi_updt && verbose
        @info @sprintf("Estimated global shrinkage parameter: %1.2e", phi_est)
    end

    return Array(beta_est)
end
