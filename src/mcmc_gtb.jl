# Ported from PRCcs/src/mcmc_gtb.py

function mcmc(a, b, phi, sst_df, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, beta_std, seed; verbose=false, profile=false, out_path)
    # seed
    if seed !== nothing
        Random.seed!(seed)
    end

    # derived stats
    beta_mrg = copy(sst_df.BETA)
    maf = copy(sst_df.MAF)
    n_pst = (n_iter-n_burnin)/thin
    p = length(sst_df.SNP)
    n_blk = length(ld_blk)

    # initialization
    beta = zeros(p)
    psi = ones(p)
    sigma = 1.0
    if phi === nothing
        phi = 1.0
        phi_updt = true
    else
        phi_updt = false
    end

    beta_est = zeros(p)
    psi_est = zeros(p)
    sigma_est = 0.0
    phi_est = 0.0

    for kk in 1:n_blk
        @assert issymmetric(ld_blk[kk])
        ld_blk[kk] = Symmetric(ld_blk[kk])
    end

    if profile
        Profile.start_timer()
    end

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
                dinvt = Symmetric(ld_blk[kk] .+ Diagonal(1.0 ./ psi[idx_blk]))
                dinvt_chol = cholesky(dinvt).U
                beta_tmp = (transpose(dinvt_chol) \ beta_mrg[idx_blk]) .+ sqrt(sigma/n) .* randn(length(idx_blk))
                beta[idx_blk] = dinvt_chol \ beta_tmp
                quad += dot(transpose(beta[idx_blk]) * dinvt, beta[idx_blk])
                mm += blk_size[kk]
            end
        end

        err = max(n/2.0*(1.0-2.0*sum(beta.*beta_mrg)+quad), n/2.0*sum(beta .^ 2 ./ psi))
        sigma = 1.0/rand(Gamma((n+p)/2.0, 1.0/err))

        delta = rand.(Gamma.(a+b, 1.0 ./ (psi .+ phi)))

        for jj in 1:p
            psi[jj] = gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]^2/sigma)
        end
        psi[psi .> 1] .= 1.0

        if phi_updt
            w = rand(Gamma(1.0, 1.0/(phi+1.0)))
            phi = rand(Gamma(p*b+0.5, 1.0/(sum(delta)+w)))
        end

        if profile && (itr == n_burnin)
            # restart profiler
            Profile.stop_timer()
            Profile.clear()
            Profile.start_timer()
        end

        # posterior
        if (itr>n_burnin) && (itr % thin == 0)
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst
        end
    end

    if profile
        Profile.stop_timer()
        open("$(out_path)_chr$chrom.cpuprofile", "w") do io
            Profile.print(io; C=true, mincount=100, maxdepth=100)
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

    return beta_est
end
