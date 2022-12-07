/*
 *  Latent Variable Model v4.0
 *    - 2021/03/01 @Yutaka Osada
 */

#include <TMB.hpp>

template <class Type>
Type penalize_mean (const Type weight, const matrix<Type> &x);

template <class Type>
Type penalize_cov (const Type weight, const matrix<Type> &x);


template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_MATRIX(Y);
    DATA_MATRIX(X);
    DATA_VECTOR(Penalty);
    
    DATA_SPARSE_MATRIX(G0);
    DATA_SPARSE_MATRIX(G1);
    DATA_SPARSE_MATRIX(G2);
    
    PARAMETER_MATRIX(beta);
    PARAMETER_VECTOR(psiSv);
    PARAMETER_VECTOR(psiRv);
    PARAMETER_VECTOR(alpha);
    PARAMETER_MATRIX(omegaS);
    PARAMETER_MATRIX(scoreR);
    PARAMETER(log_sigma);
    PARAMETER_VECTOR(log_kappa);
    
    Type nll = Type(0.0);
    const int nS = Y.cols();      // number of species
    const int nL = Y.rows();      // number of locations
    const int nE = beta.cols();   // number of environmental factors
    const int nF = omegaS.cols(); // number of spatial factors
    const int nG = scoreR.cols(); // number of local factors
    
    //* Make loading matrix (psi) as matrix upper triangle */
    matrix<Type> psiS(nS, nF);
    matrix<Type> psiR(nS, nG);
    int countS = 0, countR = 0;
    for (int i = 0; i < nS; ++i)
    {
        for (int k = 0; k < nF; ++k)
        {
            psiS(i, k) = i >= k ? psiSv(countS++): Type(0.0);
        }
        
        for (int k = 0; k < nG; ++k)
        {
            psiR(i, k) = i >= k ? psiRv(countR++): Type(0.0);
        }
    }
    
    //* latent variables (score) */
    vector<Type> kappa(nF);
    matrix<Type> scoreS(nL, nF);
    for (int k = 0; k < nF; ++k)  // Spatial latent variables
    {
        // Score Matrix - Gaussian Markov random fields
        kappa(k) = exp(log_kappa(k));
        Type kappa2 = kappa(k) * kappa(k);
        Type kappa4 = kappa2 * kappa2;
        Eigen::SparseMatrix<Type> Q = kappa4 * G0 + Type(2) * kappa2 * G1 + G2;
        nll += density::GMRF(Q)(omegaS.col(k));
        
        for (int j = 0; j < nL; ++j)
        {
            scoreS(j, k) = Type(3.545) * omegaS(j, k) * kappa(k);
            // assure var(score) == 1.0 if nu = 1.0
        }
    }
    
    for (int k = 0; k < nG; ++k)  // Local latent variables
    {
        for (int j = 0; j < nL; ++j)
        {
            nll -= dnorm(scoreR(j, k), Type(0), Type(1), true);
        }
    }
    
    Type sigma = exp(log_sigma);
    for (int j = 0; j < nL; ++j)
    {
        nll -= dnorm(alpha(j), Type(0), sigma, true);
    }
    
    //* Likelihood contribution from observations */
    for (int i = 0; i < nS; ++i)
    {
        for (int j = 0; j < nL; ++j)
        {
            Type mu = alpha(j);
            for (int k = 0; k < nE; ++k) mu += beta(i, k) * X(j, k);
            for (int k = 0; k < nF; ++k) mu += psiS(i, k) * scoreS(j, k);
            for (int k = 0; k < nG; ++k) mu += psiR(i, k) * scoreR(j, k);
            nll -= dbinom(Y(j, i), Type(1), pnorm(mu), true);
        }
    }
    
    // Penalty
    nll += penalize_mean(Penalty(0), scoreS);
    nll += penalize_mean(Penalty(1), scoreR);
    nll += penalize_cov (Penalty(2), scoreS);
    nll += penalize_cov (Penalty(3), scoreR);
    
    // Report
    REPORT(sigma);
    REPORT(kappa);
    REPORT(beta);
    REPORT(psiS);
    REPORT(psiR);
    REPORT(alpha);
    REPORT(scoreS);
    REPORT(scoreR);
    
    return nll;
}


//*----------------------------------------------------------------------*/

template <class Type>
Type penalize_mean (const Type weight, const matrix<Type> &x)
{
    Type op = Type(0);
    if (weight > Type(0))
    {
        const int K = x.cols();
        for (int k = 0; k < K; ++k)
        {
            op += weight * pow(x.col(k).mean(), 2);
        }
    }
    return op;
}


template <class Type>
Type penalize_cov (const Type weight, const matrix<Type> &x)
{
    Type op = Type(0);
    if (weight > Type(0) && x.cols() >= 2)
    {
        const int J = x.rows();
        const int K = x.cols();
        vector<Type> mean(K);
        for (int k = 0; k < K; ++k) mean(k) = x.col(k).mean();
        
        for (int k = 0; k < K; ++k)
        {
            for (int m = k + 1; m < K; ++m)
            {
                Type cov = 0; 
                for (int j = 0; j < J; ++j)
                {
                    cov += (x(j, k) - mean(k)) * (x(j, m) - mean(m));
                }
                cov /= Type(J);
                op += weight * pow(cov, 2);
            }
        }
    }
    return op;
}

//* End of file */
