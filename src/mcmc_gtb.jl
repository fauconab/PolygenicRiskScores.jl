# Ported from PRCcs/src/mcmc_gtb.py

function mcmc(a, b, phi, sst_df, n, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, seed)
    println("Initiating MCMC")

    # seed
    if seed !== nothing
        Random.seed!(seed)
    end

    # derived stats
    # FIXME: array.T
    beta_mrg = sst_df.BETA.T
    maf = sst_df.MAF.T
    n_pst = (n_iter-n_burnin)/thin
    p = length(sst_df.SNP)
    n_blk = length(ld_blk)

    # initialization
    beta = zeros(p,1)
    psi = ones(p,1)
    sigma = 1.0
    if phi === nothing
        phi = 1.0
        phi_updt = true
    else
        phi_updt = false
    end

    beta_est = zeros(p,1)
    psi_est = zeros(p,1)
    sigma_est = 0.0
    phi_est = 0.0

    # MCMC
    for itr in 1:n_iter
        if itr % 100 == 0
            println("(MCMC) iteration $itr")
        end

        mm = 0; quad = 0.0
        for kk in 1:n_blk-1
            if blk_size[kk] == 0
                continue
            else
                idx_blk = range(mm,mm+blk_size[kk])
                dinvt = ld_blk[kk]+sp.diag(1.0/psi[idx_blk].T[0])
                dinvt_chol = linalg.cholesky(dinvt)
                beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[idx_blk], trans='T') + sp.sqrt(sigma/n)*random.randn(len(idx_blk),1)
                beta[idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                quad += sp.dot(sp.dot(beta[idx_blk].T, dinvt), beta[idx_blk])
                mm += blk_size[kk]
            end
        end

        err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))
        sigma = 1.0/random.gamma((n+p)/2.0, 1.0/err)

        delta = random.gamma(a+b, 1.0/(psi+phi))

        for jj in 1:p-1
            psi[jj] = gigrnd.gigrnd(a-0.5, 2.0*delta[jj], n*beta[jj]**2/sigma)
        end
        psi[psi>1] = 1.0

        if phi_updt
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))
        end

        # posterior
        if (itr>n_burnin) && (itr % thin == 0)
            beta_est = beta_est + beta/n_pst
            psi_est = psi_est + psi/n_pst
            sigma_est = sigma_est + sigma/n_pst
            phi_est = phi_est + phi/n_pst
        end

    # convert standardized beta to per-allele beta
    if beta_std == "False"
        beta_est /= sqrt(2.0*maf*(1.0-maf))
    end

    # write posterior effect sizes
    if phi_updt
        eff_file = out_dir * "_pst_eff_a%d_b%.1f_phiauto_chr%d.txt" % (a, b, chrom)
    else
        eff_file = out_dir * "_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt" % (a, b, phi, chrom)
    end

    open(eff_file, 'w') do ff
        for snp, bp, a1, a2, beta in zip(sst_df.SNP, sst_df.BP, sst_df.A1, sst_df.A2, beta_est)
            write(ff, "%d\t%s\t%d\t%s\t%s\t%.6e\n" % (chrom, snp, bp, a1, a2, beta))
        end
    end

    # print estimated phi
    if phi_updt
        println("Estimated global shrinkage parameter: %1.2e") #FIXME: % phi_est )
    end

    println("Completed MCMC")
end
