/*  MAGEE : An R Package for Mixed Model Association Test for GEne-Environment Interaction
 *  Copyright (C) 2020  Xinyu Wang, Han Chen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2003-2004  The R Foundation
 *  Copyright (C) 1998--2013  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <bzlib.h>
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include "read_bgen.h"
#include "zstd-1.4.5/lib/zstd.h"
using namespace std;
using namespace arma;
using namespace Rcpp;


typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

#define DBL_EPSILON 2.2204460492503131e-16;


arma::vec MAF_weights_beta_fun (arma::vec freq, double beta1, double beta2) {
  
  Rcpp::Environment pkg = Environment::namespace_env("stats");
  Rcpp::Function dbeta_r = pkg["dbeta"];
  
  vector<double> out(freq.size());
  for(size_t f = 0; f < freq.size(); f++) {
    if (freq(f) > 0.5) {
      freq(f) = 1 - freq(f);
    }
    if (freq(f) <= 0) {
      freq(f) = 0;
    } else {
      Rcpp::NumericVector tmp_out = dbeta_r(Rcpp::_["x"] = freq(f), Rcpp::_["shape1"] = beta1, Rcpp::_["shape2"] = beta2);
      freq(f) = tmp_out[0];
    }
    
  }
  
  return freq;
}

double fisher_pval (arma::vec p) {

  uvec idx = find((p > 0) && (p <= 1));
  if (idx.size() == 0) {
    return NA_REAL;
  }
  
  p = p.rows(idx);
  double pval = Rf_pchisq(-2*sum(log(p)), 2*p.size(), 0, 0);
  return pval;
}

double Q_pval (double Q, arma::vec lambda, string method) {
   Rcpp::Environment pkg = Environment::namespace_env("CompQuadForm");
   Rcpp::Environment mge = Environment::namespace_env("MAGEE");
   Rcpp::Function davies_r = pkg["davies"];
   Rcpp::Function liu_r = pkg["liu"];
   Rcpp::Function pKuonen_r = mge["pKuonen"];
   double pval = 0.0;
   if(method == "davies") {
     Rcpp::List tmpList = davies_r(Rcpp::_["q"] = Q, Rcpp::_["lambda"] = lambda, Rcpp::_["acc"] = 1e-6);
     pval = Rcpp::as<double>(tmpList["Qq"]);
     double ifault = Rcpp::as<double>(tmpList["ifault"]);
     if ((ifault > 0) || (pval <= 1e-5) || (pval >= 1)) {
       method = "kuonen";
     }
   }
   if(method == "kuonen") {
     Rcpp::NumericVector tmpVec = pKuonen_r(Rcpp::_["x"] = Q, Rcpp::_["lambda"] = lambda);
     if(tmpVec == R_NilValue) {
        method = "liu";
      }
   }
   if(method == "liu") {
     if (sum(lambda) != 0) {
       Rcpp::NumericVector tmpVec = liu_r(Rcpp::_["q"] = Q, Rcpp::_["lambda"] = lambda);
       pval = tmpVec[0];
     } else {
       return NA_REAL;
     }
   }
   
   return pval;
}

double quad_pval (arma::mat U_mat, arma::mat V_mat, string method) {
   vec U2 = U_mat % U_mat;
   double Q = arma::sum(U2);
   vec lambda;
   std::ostream nullstream(0);
   arma::set_cerr_stream(nullstream);
   arma::eig_sym(lambda, V_mat);
   uvec lmb_idx = find(lambda > 0);
   lambda = lambda.rows(lmb_idx);
   double pval = Q_pval(Q, lambda, method = method);
   return(pval);
}


extern "C" 
{
  
  
  SEXP magee_bgen13(SEXP res_in, SEXP nullObj_in, SEXP bgenfile_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, 
                    SEXP isNullP_in, SEXP n_groups_in, SEXP mf_in, SEXP mv_in, SEXP if_in, SEXP iv_in, SEXP jd_in, SEXP jf_in, SEXP jv_in,
                    SEXP groupInfo_in, SEXP groupStart_in, SEXP groupEnd_in, SEXP weight_in, SEXP method_in, SEXP beta1_in, SEXP beta2_in,
                    SEXP select_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
       
    try {
      int n_groups = Rcpp::as<int>(n_groups_in);
      Rcpp::StringVector group_name(n_groups);
      Rcpp::NumericVector n_variants(n_groups);
      Rcpp::NumericVector miss_min(n_groups, 100.0);
      Rcpp::NumericVector miss_mean(n_groups);
      Rcpp::NumericVector miss_max(n_groups, -100.0);
      Rcpp::NumericVector freq_min(n_groups, 100.0);
      Rcpp::NumericVector freq_mean(n_groups);
      Rcpp::NumericVector freq_max(n_groups, -100.0);
      vec freq_strata_min = zeros(n_groups);
      vec freq_strata_max = zeros(n_groups);

      bool MF = Rcpp::as<bool>(mf_in);
      bool MV = Rcpp::as<bool>(mv_in);
      bool IF = Rcpp::as<bool>(if_in);
      bool IV = Rcpp::as<bool>(iv_in);
      bool JD = Rcpp::as<bool>(jd_in);
      bool JF = Rcpp::as<bool>(jf_in);
      bool JV = Rcpp::as<bool>(jv_in);

      Rcpp::NumericVector MF_pval;
      Rcpp::NumericVector MV_pval;
      Rcpp::NumericVector IF_pval;
      Rcpp::NumericVector IV_pval;
      Rcpp::NumericVector JD_pval;
      Rcpp::NumericVector JF_pval;
      Rcpp::NumericVector JV_pval;

      if (MV || JV) {MV_pval = NumericVector(n_groups, NumericVector::get_na());}
      if (IV || JV) {IV_pval = NumericVector(n_groups, NumericVector::get_na());}
      if (JV) {JV_pval = NumericVector(n_groups, NumericVector::get_na());}
      if (MF || JF || JD) {MF_pval = NumericVector(n_groups, NumericVector::get_na());}
      if (IF || JF || JD) {IF_pval = NumericVector(n_groups, NumericVector::get_na());}
      if (JF) {JF_pval = NumericVector(n_groups, NumericVector::get_na());}
      if (JD) {JD_pval = NumericVector(n_groups, NumericVector::get_na());}


      arma::mat P;
      arma::sp_mat Sigma_i;
      arma::sp_mat Sigma_iX;
      arma::sp_mat cov;

      Rcpp::List null_obj(nullObj_in);
      arma::mat E = as<arma::mat>(null_obj["E"]);
      bool isNullP = Rcpp::as<bool>(isNullP_in);
      if (!isNullP) {
        
        P = as<arma::mat>(null_obj["P"]);
       (void)Sigma_i; (void)Sigma_iX; (void)cov;
       } else {
         Sigma_i = as<arma::sp_mat>(null_obj["Sigma_i"]);
         Sigma_iX = as<arma::sp_mat>(null_obj["Sigma_iX"]);
         cov = as<arma::sp_mat>(null_obj["cov"]);
         (void)P;
       }

      Rcpp::NumericVector res_r(res_in);
      arma::vec res(res_r.begin(), res_r.size(), false);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const string miss_method = Rcpp::as<string>(miss_method_in);
      const double beta1 = Rcpp::as<double>(beta1_in);
      const double beta2 = Rcpp::as<double>(beta2_in);
      string method = Rcpp::as<string>(method_in);
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n);
      vector <string> biminfo;
      double gmean, geno, gmax, gmin;
      const double tol = 1e-5;
      size_t ncount, nmiss;

      Rcpp::DataFrame groupInfo(groupInfo_in);
      vector<string> groupName = Rcpp::as<vector<string>>(groupInfo["group"]);
      vector<int>groupStart = Rcpp::as<vector<int>>(groupStart_in);
      vector<int>groupEnd = Rcpp::as<vector<int>>(groupEnd_in);
      vector<long long unsigned int> bytes = Rcpp::as<vector<long long unsigned int>>(groupInfo["byte"]);
      vec weights = as<arma::vec>(weight_in);
      groupInfo.erase(0, groupInfo.size()-1);
      uint maxLA = 65536;
      std::vector<uchar> zBuf12;
      std::vector<uchar> shortBuf12;
      uint compression = Rcpp::as<uint>(compression_in);
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];

      FILE* fp = fopen(bgenfile.c_str(), "rb");

      int ret;
      for (int grp = 0; grp < n_groups; grp++) {

        int grpSize = (groupEnd[grp] - (groupStart[grp]-1));
        int npbidx = 0;
        mat G(n, grpSize);
        vec weight1 = weights.rows(groupStart[grp] - 1, groupEnd[grp]-1);
        vec freq(groupEnd[grp] - (groupStart[grp]-1));
        vec weight(groupEnd[grp] - (groupStart[grp]-1));
        group_name[grp] = groupName[groupStart[grp]-1];
        double missMean = 0.0;
        double freqMean = 0.0;
        for (int m = (groupStart[grp]-1); m < groupEnd[grp]; m++) {
            fseek(fp, bytes[m], SEEK_SET);

            uint16_t LS;
            ret = fread(&LS, 2, 1, fp);
            ret = fread(snpID, 1, LS, fp);
            snpID[LS] = '\0';

            uint16_t LR;
            ret = fread(&LR, 2, 1, fp);
            ret = fread(rsID, 1, LR, fp);
            rsID[LR] = '\0';

            uint16_t LC;
            ret = fread(&LC, 2, 1, fp);
            ret = fread(chrStr, 1, LC, fp);
            chrStr[LC] = '\0';

            uint32_t physpos;
            ret = fread(&physpos, 4, 1, fp);
            string physpos_tmp = to_string(physpos);

            uint16_t LKnum;
            ret = fread(&LKnum, 2, 1, fp);
            if (LKnum != 2) {
              Rcout << "Error reading BGEN file: There are non-bi-allelic variants. \n"; return R_NilValue;
            }

            uint32_t LA;
            ret = fread(&LA, 4, 1, fp);
            ret = fread(allele1, 1, LA, fp);
            allele1[LA] = '\0';

            uint32_t LB;
            ret = fread(&LB, 4, 1, fp);
            ret = fread(allele0, 1, LB, fp);
            allele0[LB] = '\0';

            uint cLen; ret = fread(&cLen, 4, 1, fp);

            uchar* bufAt;
            if (compression == 1) {
              zBuf12.resize(cLen - 4);
              uint dLen; ret = fread(&dLen, 4, 1, fp);
              ret = fread(&zBuf12[0], 1, cLen - 4, fp);
              shortBuf12.resize(dLen);
              uLongf destLen = dLen;

              if (uncompress(&shortBuf12[0], &destLen, &zBuf12[0], cLen-4) != Z_OK || destLen != dLen) {
                Rcout << "Error reading bgen file: Decompressing variant block failed with zLib. \n"; return R_NilValue;
              }
              bufAt = &shortBuf12[0];
            }
            else if (compression == 2) {
              zBuf12.resize(cLen - 4);
              uint dLen; ret = fread(&dLen, 4, 1, fp);
              ret = fread(&zBuf12[0], 1, cLen - 4, fp);
              shortBuf12.resize(dLen);

              uLongf destLen = dLen;
              size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
              if (ret > destLen) {
                if (ZSTD_isError(ret)) {
                  Rcout << "Error reading bgen file: Decompressing genotype block failed. \n"; return R_NilValue;
                }
              }
              bufAt = &shortBuf12[0];
            }
            else {
              zBuf12.resize(cLen);
              ret = fread(&zBuf12[0], 1, cLen, fp);
              bufAt = &zBuf12[0];
            }

            uint32_t N; memcpy(&N, bufAt, sizeof(int32_t));
            uint16_t K; memcpy(&K, &(bufAt[4]), sizeof(int16_t));
            if (K != 2) {Rcout << "Error reading bgen file: There are variants with more than 2 alleles. \n"; return R_NilValue;}
            const uint32_t min_ploidy = bufAt[6];
            if (min_ploidy != 2) {Rcout << "Error reading bgen file: Minimum ploidy should be 2. \n"; return R_NilValue;}
            const uint32_t max_ploidy = bufAt[7];
            if (max_ploidy != 2) {Rcout << "Error reading bgen file: Maximum ploidy should be 2. \n"; return R_NilValue;}

            const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
            const unsigned char* probs_start = &(bufAt[10 + N]);
            const uint32_t is_phased = probs_start[-2];
            if (is_phased < 0 || is_phased > 1) {Rcout << "Error reading bgen file: Phased value must be 0 or 1. \n"; return R_NilValue;}

            const uint32_t B = probs_start[-1];

            if (B != 8 && B!= 16 && B !=24 && B != 32) {
              Rcout << "Error reading bgen file: Bits to store probabilities must be 8, 16, 24, or 32. \n"; return R_NilValue;
            }

            const uintptr_t numer_mask = (1U << B) - 1;
            const uintptr_t probs_offset = B / 8;

            gmean=0.0;
            gmax=-100.0;
            gmin=100.0;
            nmiss=0;
            ncount = 0;

            if (!is_phased) {
              for (size_t i = 0; i < N; i++) {
                const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
                uintptr_t numer_aa;
                uintptr_t numer_ab;

                if (missing_and_ploidy == 2){
                  Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
                  probs_start += (probs_offset * 2);

                } else if (missing_and_ploidy == 130){
                  probs_start += (probs_offset * 2);
                  gmiss[select[ncount]-1] = 1;
                  nmiss++;
                  ncount++;
                  continue;

                } else {
                  Rcout << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n"; return R_NilValue;
                }


                if (select[ncount] > 0){
                  double p11 = numer_aa / double(1.0 * numer_mask);
                  double p10 = numer_ab / double(1.0 * numer_mask);
                  geno = 2 * (1 - p11 - p10) + p10;
                  gmiss[select[ncount]-1] = 0;
                  g[select[ncount]-1] = geno;
                  gmean += geno;
                  if (geno > gmax) { gmax = geno; }
                  if (geno < gmin) { gmin = geno; }
                }
                ncount++;

              }

            } else {
              for (size_t i = 0; i < N; i++) {
                const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
                uintptr_t numer_aa;
                uintptr_t numer_ab;

                if (missing_and_ploidy == 2){
                  Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
                  probs_start += (probs_offset * 2);

                } else if (missing_and_ploidy == 130){
                  probs_start += (probs_offset * 2);
                  gmiss[select[ncount]-1] = 1;
                  nmiss++;
                  ncount++;
                  continue;

                } else {
                  Rcout << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n"; return R_NilValue;
                }

                if (select[ncount] > 0){
                  double p11 = numer_aa / double(1.0 * numer_mask);
                  double p10 = numer_ab / double(1.0 * numer_mask);
                  geno = double(1.0 * missing_and_ploidy) - (p11 + p10);
                  gmiss[select[ncount]-1] = 0;
                  g[select[ncount]-1] = geno;
                  gmean += geno;
                  if (geno > gmax) { gmax = geno; }
                  if (geno < gmin) { gmin = geno; }
                }
                ncount++;

              }
            }
            gmean/=(double)(n-nmiss);
            if (nmiss > 0) {
              for (size_t j=0; j<n; ++j) {
                if (gmiss[j]==1) {
                  if (miss_method == "impute2mean") {
                    g[j] = gmean;
                  } else {
                    g[j] = 0.0;
                  }
                }
              }
            }

            gmean /= 2.0; // convert mean to allele freq
            if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
              continue;

            } else {
              G.col(npbidx) = g;
              double miss = (double)(nmiss/n);
              freqMean += gmean;
              missMean += miss;
              freq(npbidx) = gmean;
              miss_min[grp]  = (miss < miss_min[grp]) ? miss : miss_min[grp];
              miss_max[grp]  = (miss > miss_max[grp]) ? miss : miss_max[grp];
              freq_min[grp] = (gmean < freq_min[grp]) ? gmean : freq_min[grp];
              freq_max[grp] = (gmean > freq_max[grp]) ? gmean : freq_max[grp];
              weight(npbidx) = weight1(npbidx);
              npbidx++;
            }


        }
        if (npbidx == 0) {
          miss_min[grp] = NA_REAL;
          miss_mean[grp] = NA_REAL;
          miss_max[grp] = NA_REAL;
          freq_min[grp] = NA_REAL;
          freq_mean[grp] = NA_REAL;
          freq_max[grp] = NA_REAL;
          continue;
        }
        n_variants[grp] = npbidx;
        miss_mean[grp] = missMean / npbidx;
        freq_mean[grp] = freqMean / npbidx;
        G.reshape(n, npbidx);
        freq.reshape(npbidx, 1);
        weight.reshape(npbidx, 1);
        arma::vec U = G.t() * res;
        sp_mat PG;
        sp_mat Gsp(G);
        if (!isNullP) {
          PG = P.t() * G;
        } else {
          sp_mat GSigma_iX =  Gsp.t() * Sigma_iX;
          PG = (Sigma_i.t() * Gsp) - (Sigma_iX * (GSigma_iX * cov.t()).t());
        }
        mat V = G.t() * PG;
        mat K;
        mat IV_U;
        mat IV_V;
        if(IV || IF || JV || JF || JD) {
          for (size_t ei_idx = 0; ei_idx < E.n_cols; ei_idx++) {
            mat tmpK = G.each_col() % E.col(ei_idx);
            K = join_horiz(K, tmpK);
          }
          mat KPK;
          if (!isNullP) {
            KPK = K.t() * (P.t() * K);
          } else {
            mat KSigma_iX = K.t() * Sigma_iX;
            KPK = (K.t() * (Sigma_i.t() * K)) - (KSigma_iX * (KSigma_iX * cov.t()).t());
          }
          mat V_i = inv(V);
          mat KPG = K.t() * PG;
          IV_U = (K.t() * res) - ((KPG * V_i.t()) * U);
          IV_V = KPK - ((KPG * V_i.t()) * KPG.t());
        }
        vec tmpWeight = MAF_weights_beta_fun(freq, beta1, beta2);
        weight = weight % tmpWeight;
        U = U % weight;
        V = (V.each_col() % weight).t();
        V = V.each_col() % weight;
        if(IV || IF || JV || JF || JD) {
          vec tmpWeight2;
          if (E.n_cols > 1) {
            for (size_t tmp = 0; tmp < (E.n_cols); tmp++){
              tmpWeight2 = join_vert(tmpWeight2, weight);
            }
            IV_U = IV_U.each_col() % tmpWeight2;
            IV_V = (IV_V.each_col() % tmpWeight2).t();
            IV_V = IV_V.each_col() % tmpWeight2;
          } else {
            IV_U = IV_U.each_col() % weight;
            IV_V = (IV_V.each_col() % weight).t();
            IV_V = IV_V.each_col() % weight;
          }
        }

        if (MV || JV) { MV_pval[grp] = quad_pval(U,  V, method);}
        if (IV || JV) { IV_pval[grp] = quad_pval(IV_U,  IV_V, method);}
        vec MV_IV_pval(2); MV_IV_pval(0) = MV_pval[grp]; MV_IV_pval(1) = IV_pval[grp];
        if (JV) { JV_pval[grp] = fisher_pval(MV_IV_pval);}
        double MF_Bp = 0.0;
        double MF_p = 0.0;
        if (MF || JF || JD) {
          vec MF_BU(npbidx); MF_BU.fill(sum(U));
          vec MF_BV(npbidx); MF_BV.fill(accu(V));

          MF_Bp = Rf_pchisq((MF_BU(0) * MF_BU(0))/MF_BV(0), 1, 0, 0);
          vec V_rowSums = sum(V, 1);

          mat MF_U = U - ((V_rowSums % MF_BU) / MF_BV);
          mat MF_V = V_rowSums * V_rowSums.t();
          MF_V = V - (MF_V.each_col() / MF_BV);

          double eps = DBL_EPSILON;
          if ((accu(abs(MF_V)) / (npbidx * npbidx)) < sqrt(eps)) {
            MF_p = NA_REAL;
          } else {
            MF_p = quad_pval(MF_U, MF_V, method);
          }
          vec MFBp_MFp_pval(2); MFBp_MFp_pval(0) = MF_Bp; MFBp_MFp_pval(1) = MF_p;
          MF_pval[grp] = fisher_pval(MFBp_MFp_pval);
        }
        double IF_Bp = 0.0;
        double IF_p = 0.0;
        if (IF || JF || JD) {
          vec IF_BU(npbidx * E.n_cols); IF_BU.fill(accu(IV_U));
          vec IF_BV(npbidx * E.n_cols); IF_BV.fill(accu(IV_V));

          IF_Bp = Rf_pchisq((IF_BU(0) * IF_BU(0))/IF_BV(0), 1, 0, 0);
          vec IV_rowSums = sum(IV_V, 1);
          mat IF_U = IV_U - ((IV_rowSums % IF_BU) / IF_BV);
          mat IF_V = IV_rowSums * IV_rowSums.t();
          IF_V = IV_V - (IF_V.each_col() / IF_BV);

          double eps = DBL_EPSILON;
          if ((accu(abs(IF_V)) / (npbidx * npbidx)) < sqrt(eps)) {
            IF_p = -9.0;
          } else {
            IF_p = quad_pval(IF_U, IF_V, method);
          }
          vec IFBp_IFp_pval(2); IFBp_IFp_pval(0) = IF_Bp; IFBp_IFp_pval(1) = IF_p;
          IF_pval[grp] = fisher_pval(IFBp_IFp_pval);
        }

        vec jf_tmp(4); jf_tmp(0) = MF_Bp; jf_tmp(1) = MF_p; jf_tmp(2) = IF_Bp; jf_tmp(3) = IF_p;
        if (JF) {JF_pval[grp] = fisher_pval(jf_tmp);}
        vec jd_tmp(2); jd_tmp(0) = MF_pval[grp]; jd_tmp(1) = IF_pval[grp];
        if (JD) {JD_pval[grp] = fisher_pval(jd_tmp);}

      }


      (void)ret;
      delete [] snpID;
      delete [] rsID;
      delete [] chrStr;
      delete [] allele0;
      delete [] allele1;
      fclose(fp);
      return(Rcpp::DataFrame::create(Named("group") = group_name,
                                     Named("n.variants") = n_variants,
                                     Named("miss.min") = miss_min,
                                     Named("miss.mean") = miss_mean,
                                     Named("miss.max") = miss_max,
                                     Named("freq.min") = freq_min,
                                     Named("freq.mean") = freq_mean,
                                     Named("freq.max") = freq_max,
                                     Named("MV.pval") = MV_pval,
                                     Named("MF.pval") = MF_pval,
                                     Named("IV.pval") = IV_pval,
                                     Named("IF.pval") = IF_pval,
                                     Named("JV.pval") = JV_pval,
                                     Named("JF.pval") = JF_pval,
                                     Named("JD.pval") = JD_pval));
    } 
    catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;

  }  
  
} // end of extern "C"
