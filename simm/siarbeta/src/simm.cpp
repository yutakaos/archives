/*------------------------------------------------------------------------------------------#
 * Beta-dependent Stable Isotope Mixing Model
 *------------------------------------------------------------------------------------------#
 */

#ifndef _simm_cpp_
#define _simm_cpp_

/* Header(s) */
#include <Rcpp.h>
#include <vector>
#include <utils.hpp>

typedef double num_t;

class Simm
{
    
private:
    
    std::vector<std::vector<num_t>> DataMIX;
    std::vector<num_t> alpha;
    num_t tempL, tempP, lp_const;
    
    int nR, nE;
    MIXTURE<num_t> mix;
    SOURCE <num_t> src;
    
    int  ERROR_TYPE;
    bool SET_MX, SET_SX, SET_SQ, SET_SD, SET_PR;
    
public:
    
    Simm() :
        tempL(1), tempP(1), ERROR_TYPE(3),
        SET_MX(false), SET_SX(false), SET_SQ(false), SET_SD(false), SET_PR(false)
    {};
    
    void set_srcX_summary (const Rcpp::NumericMatrix &srcX)
    {
        nR = srcX.nrow();
        nE = srcX.ncol() / 2;
        utils::set_summary_data(srcX, &src.Xm, &src.Xv);
        SET_SX = true;
    }
    
    void set_srcQ_summary (const Rcpp::NumericMatrix &srcQ)
    {
        if (!SET_SX) Rcpp::stop("Unspecified source data.");
        if (srcQ.nrow() !=     nR) Rcpp::stop("Different number of sources.");
        if (srcQ.ncol() != 2 * nE) Rcpp::stop("Different number of elements.");
        
        std::vector<std::vector<num_t>> dummy;
        utils::set_summary_data(srcQ, &src.Qm, &dummy);
        SET_SQ = true;
    }
    
    void set_srcD_summary (const Rcpp::NumericMatrix &srcD)
    {
        if (!SET_SX) Rcpp::stop("Unspecified source data.");
        if (srcD.nrow() !=     nR) Rcpp::stop("Different number of sources.");
        if (srcD.ncol() != 2 * nE) Rcpp::stop("Different number of elements.");
        
        utils::set_summary_data(srcD, &src.Dm, &src.Dv);
        SET_SD = true;
    }
    
    void set_mixX (const Rcpp::NumericMatrix &mixX)
    {
        if (!SET_SX) Rcpp::stop("Unspecified source data.");
        if (mixX.ncol() != nE) Rcpp::stop("Different number of elements.");
        
        std::vector<std::vector<num_t>>().swap(DataMIX);
        int nD = mixX.nrow();
        DataMIX.resize(nD, std::vector<num_t>(nE));
        for (int i = 0; i < nD; ++i)
        {
            for (int j = 0; j < nE; ++j)
            {
                DataMIX[i][j] = mixX(i, j);
            }
        }
        SET_MX = true;
    }
    
    void set_priors (const Rcpp::NumericVector &alpha)
    {
        if (!SET_SX) Rcpp::stop("Unspecified source data.");
        if (alpha.size() != nR) Rcpp::stop("Different number of sources.");
        
        std::vector<num_t>().swap(this->alpha);
        this->alpha.resize(nR);
        for (int i = 0; i < nR; ++i) this->alpha[i] = alpha(i);
        SET_PR = true;
    }
    
    void set_temp (const num_t tempL, const num_t tempP)
    {
        this->tempL = tempL;
        this->tempP = tempP;
    }
    
    void initialize (const int error_type)
    {
        if (!SET_MX) Rcpp::stop("Unspecified mixture data.");
        if (!SET_SQ) Rcpp::stop("Unspecified source correction data.");
        if (!SET_SD) Rcpp::stop("Unspecified source concdep data.");
        if (!SET_PR) Rcpp::stop("Unspecified alpha priors.");
        
        mix.LL = 0;
        std::vector<num_t>().swap(mix.ps);
        std::vector<num_t>().swap(mix.Xm);
        std::vector<num_t>().swap(mix.Xv);
        std::vector<num_t>().swap(mix.Rv);
        mix.ps.resize(nR);
        mix.Rv.resize(nE);
        mix.ep.resize(nE);
        initialize_error_structure(error_type);
        
        num_t sum_alpha = utils::sum(alpha);
        lp_const = std::lgamma(sum_alpha);
        for (int i = 0; i < nR; ++i)
        {
            mix.ps[i] = alpha[i] / sum_alpha;
            lp_const -= std::lgamma(alpha[i]);
        }
        utils::get_XmXv(mix.ps, src, &mix.Xm, &mix.Xv);
    }
    
    void update (const int iters)
    {
        for (int i = 0; i < iters; ++i)
        {
            update_Rs();
            update_ps(200);
        }
    }
    
    Rcpp::List output (const bool posterior)
    {
        std::vector<num_t> Rsd = mix.Rv;
        for (auto &x : Rsd) x = std::sqrt(x);
        
        if (posterior)
        {
            num_t lp = mix.LL + lp_const + PDF::dirichlet(mix.ps, alpha);
            return Rcpp::List::create(
                Rcpp::Named("ps") = Rcpp::wrap(mix.ps),
                Rcpp::Named("sd") = Rcpp::wrap(Rsd),
                Rcpp::Named("ep") = Rcpp::wrap(mix.ep),
                Rcpp::Named("lp") = lp
            );
        }
        return Rcpp::List::create(
            Rcpp::Named("ps") = Rcpp::wrap(mix.ps),
            Rcpp::Named("sd") = Rcpp::wrap(Rsd),
            Rcpp::Named("ep") = Rcpp::wrap(mix.ep)
        );
    }
    
private:
    
    void initialize_error_structure (const int error_type)
    {
        if (error_type < 0 || 3 < error_type) Rcpp::stop("Invalid error strucure.");
        if (DataMIX.size() == 1 && error_type != 0)
            Rcpp::stop("Invalid error structure for single mixture data.");
        
        ERROR_TYPE = error_type;
        if (ERROR_TYPE == 1)
        {
            for (int e = 0; e < nE; ++e)
            {
                for (int i = 0; i < nR; ++i)
                {
                    src.Xv[i][e] = 0.0;
                    src.Dv[i][e] = 0.0;
                }
            }
        }
        num_t init_Rv = (ERROR_TYPE == 1 || ERROR_TYPE == 2) ? 100 : 0;
        for (int e = 0; e < nE; ++e)
        {
            mix.Rv[e] = init_Rv;
            mix.ep[e] = 1.0;
        }
    }
    
    num_t compute_LL (
        const std::vector<num_t> &ps,
        std::vector<num_t> *Xm, std::vector<num_t> *Xv)
    {
        num_t LL = 0;
        utils::get_XmXv(ps, src, Xm, Xv);
        for (int e = 0; e < nE; ++e)
        {
            num_t sd = std::sqrt((*Xv)[e] * mix.ep[e] + mix.Rv[e]);
            for (auto &x : DataMIX)
            {
                LL += PDF::normal(x[e], (*Xm)[e], sd);
            }
        }
        return LL;
    }
    
    void update_ps (const int NP)
    {
        for (int R1 = 0; R1 < nR; ++R1) 
        {
            int R2 = int(num_t(nR - 1) * R::runif(0,1));
            if (R1 <= R2) ++R2;
            
            std::vector<num_t> Xm, Xv, ps = mix.ps;
            std::vector<num_t> lw(NP), x(NP);
            for (int i = 0; i < NP; ++i)
            {
                x[i] = R::runif(0,1);
                utils::set_ps(&ps, R1, R2, x[i]);
                num_t LP = PDF::dirichlet(ps, alpha);
                lw[i] = tempL * compute_LL(ps, &Xm, &Xv) + tempP * LP;
            }
            utils::set_ps(&mix.ps, R1, R2, x[RNG::sample_lw(lw)]);
            mix.LL = compute_LL(mix.ps, &Xm, &Xv);
            mix.Xm = Xm;
            mix.Xv = Xv;
        }
    }
    
    void update_Rs()
    {
        if (ERROR_TYPE == 0) return; // mix.Rv[e] = 0;
        
        num_t shape = 0.5 * tempL * num_t(DataMIX.size());
        std::vector<num_t> scale(nE, 0);
        for (int e = 0; e < nE; ++e)
        {
            for (auto &x : DataMIX) scale[e] += (x[e] - mix.Xm[e]) * (x[e] - mix.Xm[e]);
            scale[e] *= 0.5 * tempL;
        }
        
        if (ERROR_TYPE == 3)
        {
            for (int e = 0; e < nE; ++e)
            {
                num_t epp = RNG::igammaLT(shape, scale[e], num_t(0), mix.Xv[e] * num_t(20));
                mix.ep[e] = epp / mix.Xv[e];
            }
        }
        else
        {
            for (int e = 0; e < nE; ++e)
            {
                num_t Rvp = RNG::igammaLT(shape, scale[e], mix.Xv[e], mix.Xv[e] + num_t(100));
                mix.Rv[e] = Rvp - mix.Xv[e];
            }
        }
    }
};


RCPP_MODULE(simm)
{
    Rcpp::class_<Simm>("SIMM")
    .constructor()
    .method("set_srcX_summary", &Simm::set_srcX_summary)
    .method("set_srcQ_summary", &Simm::set_srcQ_summary)
    .method("set_srcD_summary", &Simm::set_srcD_summary)
    .method("set_mixX"  , &Simm::set_mixX)
    .method("set_priors", &Simm::set_priors)
    .method("set_temp"  , &Simm::set_temp)
    .method("init"      , &Simm::initialize)
    .method("update"    , &Simm::update)
    .method("output"    , &Simm::output);
}

#endif