# Ported from PRCcs/src/mcmc_gtb.py

function mcmc(; a, b, phi, snp_df, beta_vecs, frq_vecs, idx_vecs, sst_df, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, beta_std, meta, seed, verbose=false, profile=false)
    #TODO: fix single-ancestry mode where snp_df, beta_vecs, frq_vecs, idx_vecs are set to nothing
    # seed
    if seed !== nothing
        Random.seed!(seed)
    end

    # derived stats
    n_pop = length(n)
    n_pst = (n_iter-n_burnin)/thin
    p_tot = length(snp_df.SNP)

    p = [length(beta_vecs[pp]) for pp in 1:n_pop]
    n_blk = [length(ld_blk[pp]) for pp in 1:n_pop]
    het = [sqrt.(2.0 .* frq_vecs[pp] .* (1.0 .- frq_vecs[pp])) for pp in 1:n_pop]

    n_grp = zeros(p_tot)
    for jj in 1:p_tot
        for pp in 1:n_pop
            if jj in idx_vecs[pp]
                n_grp[jj] += 1
            end
        end
    end

    # initialization
    beta = [zeros(p[pp]) for pp in 1:n_pop]
    sigma = ones(n_pop)
    psi = ones(p_tot)
    if phi === nothing
        phi = 1.0
        phi_updt = true
    else
        phi_updt = false
    end

    beta_est = [zeros(p[pp]) for pp in 1:n_pop]
    beta_sq_est = [zeros(p[pp]) for pp in 1:n_pop]
    psi_est = zeros(p_tot)
    sigma_est = zeros(n_pop)
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

        for pp in 1:n_pop
            mm = 1; quad = 0.0
            psi_pp = psi[idx_vecs[pp]]
            for kk in 1:n_blk[pp]
                if blk_size[pp][kk] == 0
                    continue
                else
                    idx_blk = mm:(mm+blk_size[pp][kk]-1)
                    dinvt = Symmetric(ld_blk[pp][kk] .+ Diagonal(1.0 ./ psi_pp[idx_blk]))
                    dinvt_chol = cholesky(dinvt).U
                    beta_tmp = (transpose(dinvt_chol) \ beta_vecs[pp][idx_blk]) .+ sqrt(sigma[pp]/n[pp]) .* randn(length(idx_blk))
                    beta[pp][idx_blk] = dinvt_chol \ beta_tmp
                    quad += dot(transpose(beta[pp][idx_blk]) * dinvt, beta[pp][idx_blk])
                    mm += blk_size[pp][kk]
                end
            end

            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp].*beta_vecs[pp])+quad), n[pp]/2.0*sum(beta[pp] .^ 2 ./ psi_pp))
            sigma[pp] = 1.0/rand(Gamma((n[pp]+p[pp])/2.0, 1.0/err))
        end

        delta = rand.(Gamma.(a+b, 1.0 ./ (psi .+ phi)))

        xx = zeros(p_tot)
        for pp in 1:n_pop
            xx[idx_vecs[pp]] .+= n[pp] .* beta[pp] .^ 2 ./ sigma[pp]
        end

        for jj in 1:p_tot
            while true
                try
                    psi[jj] = gigrnd(a-0.5*n_grp[jj], 2.0*delta[jj], xx[jj])
                    break
                catch
                    continue
                end
            end
        end

        psi[psi .> 1] .= 1.0

        if phi_updt
            w = rand(Gamma(1.0, 1.0/(phi+1.0)))
            phi = rand(Gamma(p_tot*b+0.5, 1.0/(sum(delta)+w)))
        end

        if profile && (itr == n_burnin)
            # restart profiler
            Profile.stop_timer()
            Profile.clear()
            Profile.start_timer()
        end

        # posterior
        if (itr>n_burnin) && (itr % thin == 0)
            for pp in 1:n_pop
                beta_est[pp] = beta_est[pp] + beta[pp] ./ n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp] .^ 2 ./ n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp] ./ n_pst
            end
            phi_est = phi_est + phi/n_pst
            psi_est = psi_est + psi/n_pst
        end
    end

    if profile
        Profile.stop_timer()
    end

    # convert standardized beta to per-allele beta
    if !beta_std
        for pp in 1:n_pop
            beta_est[pp] ./= het[pp]
            beta_sq_est[pp] ./= het[pp] .^ 2
        end
    end

    if meta
        vv = zeros(p_tot)
        zz = zeros(p_tot)
        for pp in 1:n_pop
            vv[idx_vecs[pp]] .+= 1.0 ./ (beta_sq_est[pp] .- beta_est[pp] .^ 2)
            zz[idx_vecs[pp]] .+= 1.0 ./ (beta_sq_est[pp] .- beta_est[pp] .^ 2) .* beta_est[pp]
        end
        mu = zz ./ vv
    else
        mu = nothing
    end

    return beta_est, (;phi_est=phi_est, mu=mu)
end
