#ifndef SAIR_UTILS_H
#define SAIR_UTILS_H

#include <Rcpp.h>
#include <random>
#include <vector>


//* data_type */
template<typename num_t>
struct MIXTURE
{
    num_t LL;
    std::vector<num_t> ps;  // dietary proportion [R]
    std::vector<num_t> Xm;  // isotope mean       [E]
    std::vector<num_t> Xv;  // isotope variance   [E]
    std::vector<num_t> Rv;  // residual error     [E]
    std::vector<num_t> ep;  // epsilon error      [E]
};

template<typename num_t>
struct SOURCE
{
    std::vector<std::vector<num_t>> Xm; // isotope mean       [R,E]
    std::vector<std::vector<num_t>> Xv; // isotope variance   [R,E]
    std::vector<std::vector<num_t>> Qm; // concentration mean [R,E]
    std::vector<std::vector<num_t>> Dm; // TEF mean           [R,E]
    std::vector<std::vector<num_t>> Dv; // TEF variance       [R,E]
};


//* helper functions */
namespace utils
{
    template <typename num_t>
    void set_summary_data (
        const Rcpp::NumericMatrix &data,
        std::vector<std::vector<num_t>> *M, // mean data
        std::vector<std::vector<num_t>> *V) // variance data
    {
        size_t nR = data.nrow();
        size_t nE = data.ncol() / 2;
        std::vector<std::vector<num_t>>().swap(*M);
        std::vector<std::vector<num_t>>().swap(*V);
        (*M).resize(nR, std::vector<num_t>(nE));
        (*V).resize(nR, std::vector<num_t>(nE));
        for (size_t i = 0; i < nR; ++i)
        {
            for (size_t j = 0; j < nE; ++j)
            {
                size_t j2 = 2 * j;
                (*M)[i][j] = data(i, j2);
                (*V)[i][j] = data(i, j2 + 1) * data(i, j2 + 1);
            }
        }
    }
    
    template <typename num_t>
    void get_XmXv (
        const std::vector<num_t> &ps, const SOURCE<num_t> &src,
        std::vector<num_t> *Xm, std::vector<num_t> *Xv)
    {
        size_t nR = src.Xm.size();
        size_t nE = src.Xm[0].size();
        (*Xm).resize(nE);
        (*Xv).resize(nE);
        num_t xm, xv, sumpq, pq;
        for (size_t e = 0; e < nE; ++e)
        {
            xm = xv = sumpq = 0.0;
            for (size_t i = 0; i < nR; ++i)
            {
                pq = ps[i] * src.Qm[i][e];
                xm += (src.Xm[i][e] + src.Dm[i][e]) * pq;
                xv += (src.Xv[i][e] + src.Dv[i][e]) * pq * pq;
                sumpq += pq;
            } 
            (*Xm)[e] = xm / sumpq;	
            (*Xv)[e] = xv / sumpq / sumpq;
        }
    }
    
    template <typename num_t>
    num_t sum (const std::vector<num_t> &vec)
    {
        num_t sum = 0;
        for (auto x : vec) sum += x;
        return sum;
    }
    
    template <typename num_t>
    void set_ps (
        std::vector<num_t> *ps,
        const int R1, const int R2, const num_t x)
    {
        num_t sump = (*ps)[R1] + (*ps)[R2];
        (*ps)[R1] = sump * x;
        (*ps)[R2] = sump - (*ps)[R1];
    }
}


//* probability functions */
namespace RNG
{
    template <typename num_t>
    int sample_lw (const std::vector<num_t> &log_weight)
    {
        std::uniform_real_distribution<num_t> runif;
        
        int np = log_weight.size();
        num_t max_lw = log_weight[0];
        for (int i = 1; i < np; ++i)
        {
            if (max_lw < log_weight[i]) max_lw = log_weight[i];
        }
        std::vector<num_t> w(np);
        num_t sum_w = 0;
        for (int i = 0; i < np; ++i)
        {
            w[i] = exp(log_weight[i] - max_lw);
            sum_w += w[i];
        }
        num_t rand = R::runif(0, sum_w);
        for(int i = 0; i < np ; ++i)
        {
            rand -= w[i];
            if (rand < 0) return i;
        }
        return np - 1;
    }
    
    template <typename num_t>
    num_t igammaLT (
        const num_t shape, const num_t scale, const num_t min, const num_t max)
    {
        num_t logl = R::pgamma(1/max, shape, 1/scale, 1, 1);
        if (!finite(logl)) return min;
        num_t logr = R::pgamma(1/min, shape, 1/scale, 1, 1);
        num_t logu = log(R::runif(0,1)) + log(1 - exp(logl - logr)) + logr;
        num_t logm = logl < logu ? logu : logl;
        logu = log(exp(logl - logm) + exp(logu - logm)) + logm;
        return 1/R::qgamma(logu, shape, 1/scale, 1, 1);
    }
}

namespace PDF
{
    template <typename num_t>
    num_t normal (const num_t x, const num_t mean, const num_t sd)
    {
        return R::dnorm4(x, mean, sd, true); // log-scale
    }
    
    template <typename num_t>
    num_t dirichlet (const std::vector<num_t> &ps, const std::vector<num_t> &alpha)
    {
        num_t lp = 0;
        for (size_t i = 0; i < ps.size(); ++i)
        {
            lp -= (alpha[i] - 1) * std::log(ps[i]);
        }
        return lp; // log-scale
    }
}

#endif