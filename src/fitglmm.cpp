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


extern "C" 
{
  
    
  SEXP glmm_gei_bgen13(SEXP res_in, SEXP nullObj_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, 
                       SEXP ei_in, SEXP qi_in, SEXP isNullP_in, SEXP isNullEC_in,
                       SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
    try{
      
      bool isNullP = Rcpp::as<bool>(isNullP_in);
      bool isNullEC = Rcpp::as<bool>(isNullEC_in);
      int ei = Rcpp::as<int>(ei_in);
      int qi = Rcpp::as<int>(qi_in);
      
      arma::mat P;
      arma::sp_mat Sigma_i;
      arma::sp_mat Sigma_iX;
      arma::sp_mat cov;
      
      Rcpp::List null_obj(nullObj_in);
      
      arma::mat E = as<arma::mat>(null_obj["E"]);
      arma::mat EC;
      if (!isNullEC) {
        EC = as<arma::mat>(null_obj["EC"]);
      } else {
        (void)EC;
      }
      
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
      const char center = Rcpp::as<char>(center_in);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const char miss_method = Rcpp::as<char>(miss_method_in);
      
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      string outfile  = Rcpp::as<string>(outfile_in);
      
      const size_t npb = Rcpp::as<size_t>(nperbatch_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n), snp_skip = zeros<uvec>(npb);
      mat G(n, npb);
      string* tmpout = new string[npb];
      vector <string> biminfo;
      double gmean, geno, gmax, gmin;
      const double tol = 1e-5;
      size_t ncount, nmiss, npbidx = 0;
      double compute_time = 0.0;
      ofstream writefile(outfile.c_str(), ofstream::out);
      
      

      uint maxLA = 65536;
      std::vector<uchar> zBuf12;
      std::vector<uchar> shortBuf12;
      uint compression = Rcpp::as<uint>(compression_in);
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];
      
      uint begin = Rcpp::as<uint>(begin_in);
      uint end   = Rcpp::as<uint>(end_in);
      long long unsigned int byte = Rcpp::as<long long unsigned int>(pos_in);
      FILE* fp = fopen(bgenfile.c_str(), "rb");
      fseek(fp, byte, SEEK_SET);
      
      bool isMultiThread = Rcpp::as<bool>(isMultiThread_in);
      if (!isMultiThread){
        writefile << "SNP\tRSID\tCHR\tPOS\tREF\tALT\tN\tMISSRATE\tAF\tBETA.MAIN\tSE.MAIN\tPVAL.MAIN\tSTAT.INT\tPVAL.INT\tPVAL.JOINT\n";
      }
      
      int ret;
      for (uint m = begin; m < end; m++) {
        stringstream writeout;
        
        ushort LS; 
        ret = fread(&LS, 2, 1, fp);
        ret = fread(snpID, 1, LS, fp); 
        snpID[LS] = '\0';
        
        ushort LR; 
        ret = fread(&LR, 2, 1, fp);
        ret = fread(rsID, 1, LR, fp);  
        rsID[LR] = '\0';
        
        ushort LC; 
        ret = fread(&LC, 2, 1, fp);
        ret = fread(chrStr, 1, LC, fp); 
        chrStr[LC] = '\0';
        
        uint physpos; 
        ret = fread(&physpos, 4, 1, fp);
	      string physpos_tmp = to_string(physpos);
	      
        ushort LKnum; 
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
       for (size_t j=0; j<n; ++j) {
        if (gmiss[j]==1) {
          g[j] = gmean;
          if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
        }
        if (center=='c') {
          g[j] -= gmean;
        }
       }

       gmean /= 2.0; // convert mean to allele freq
       writeout << snpID << "\t" << rsID << "\t" << chrStr << "\t" << physpos_tmp << "\t" << allele1 << "\t" << allele0 << "\t" << (n-nmiss) << "\t" << (n-nmiss)/n << "\t"<< gmean << "\t";
       tmpout[npbidx] = writeout.str();
       writeout.clear();
       if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
        snp_skip[npbidx] = 1;
       } else {
        G.col(npbidx) = g;
       }
  
       npbidx++;
       
       vec V_i;
       mat IV_U;
       mat IV_V_i;
       vec BETA_MAIN;
       vec STAT_MAIN;
       vec STAT_INT;
       vec SE_MAIN;
       if((m+1 == end) || (npbidx == npb)) {
          
        
         uvec snp_idx = find(snp_skip == 0);
         if (npbidx != npb) {
           G.reshape(n, npbidx);
           snp_idx = snp_idx.rows(0,npbidx-1);
         }
         G = G.cols(snp_idx);
         
         if (G.n_cols != 0) {
           int ng = G.n_cols;
           mat K;
           for (int ei_idx = 0; ei_idx < ei; ei_idx++) {
              mat tmpK = G.each_col() % E.col(ei_idx);
              K = join_horiz(K, tmpK);
           } 
            if (!isNullEC) {
              uvec g_idx = linspace<uvec>(0, ng-1, ng);
              for (int qi_idx = 0; qi_idx < qi; qi_idx++) {
                mat tmpG = G.cols(g_idx);
                tmpG = tmpG.each_col() % EC.col(qi_idx);
                G = join_horiz(G, tmpG);
              }
            }
            
            vec U = G.t() * res;
            sp_mat Gsp(G);
            sp_mat PG;
            if (!isNullP) {
              PG = P.t() * G;
            } else {
              sp_mat GSigma_iX =  Gsp.t() * Sigma_iX;
              PG = (Sigma_i.t() * Gsp) - (Sigma_iX * (GSigma_iX * cov.t()).t());
            }
            mat GPG = (G.t() * PG)  % kron(ones(1+qi, 1+qi), mat(ng, ng, fill::eye));
            mat GPG_i = inv(GPG);
            
            if (isNullEC) {
              V_i = diagvec(GPG_i);
            } else {
              V_i = diagvec(GPG);
              V_i = V_i.rows(0, ng-1);
              double eps = DBL_EPSILON;
              for (size_t tmp_v = 0; tmp_v < V_i.size(); tmp_v++) {
                if (V_i[tmp_v] < eps) {
                  V_i[tmp_v] = 0;
                } else {
                  V_i[tmp_v] =  1 / V_i[tmp_v];
                }
              }
            }
            
            BETA_MAIN = V_i % U.rows(0,ng-1);
            SE_MAIN = sqrt(V_i);
            STAT_MAIN = BETA_MAIN % U.rows(0,ng-1);
            mat KPK;
            
            if (!isNullP) {
              KPK = K.t() * (P.t() * K);
            } else {
              mat KSigma_iX = K.t() * Sigma_iX;
              KPK = (K.t() * (Sigma_i.t() * K)) - (KSigma_iX * (KSigma_iX * cov.t()).t());
            }
            
            KPK = KPK % kron(ones(ei, ei), mat(ng, ng, fill::eye));
            mat KPG = (K.t() * PG) % kron(ones(ei, 1+qi), mat(ng, ng, fill::eye));
            if (isNullEC) {
              mat IV_U_kron = kron(ones(ei,1), mat(ng, ng, fill::eye));
              IV_U = IV_U_kron.each_col() % ((K.t() * res) - (KPG * BETA_MAIN));
              IV_V_i = inv(KPK - (KPG * (KPG * GPG_i.t()).t()));
              STAT_INT = diagvec(IV_U.t() * (IV_V_i.t() * IV_U));
              
            } else {
              mat IV_U_kron = kron(ones(ei,1), mat(ng, ng, fill::eye));
              IV_U = IV_U_kron.each_col() % ((K.t() * res) - (KPG * (GPG_i.t() * U)));
              IV_V_i = inv(KPK - (KPG * (KPG * GPG_i.t()).t()));
              STAT_INT = diagvec(IV_U.t() * (IV_V_i.t() * IV_U));
            }
          
         }
               
         int p_idx = 0;
         for(size_t j=0; j<npbidx; ++j) {
           if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
             writefile << tmpout[j] << "NA\tNA\tNA\tNA\tNA\tNA\n";
           } 
           else {
            if (V_i[j] > 0) {
              double p_m = Rf_pchisq(STAT_MAIN[p_idx], 1, 0, 0);
              double p_i = Rf_pchisq(STAT_INT[p_idx], ei, 0, 0);
              double p_j = Rf_pchisq(STAT_MAIN[p_idx] + STAT_INT[p_idx], 1+ei, 0, 0);
              writefile << tmpout[j] << BETA_MAIN[p_idx] << "\t" << SE_MAIN[p_idx] << "\t" << p_m << "\t" << STAT_INT[p_idx] << "\t" << p_i << "\t" << p_j << "\n";
            } 
            else {
              double p_i = Rf_pchisq(STAT_INT[j], ei, 0, 0);
              writefile << tmpout[j] << BETA_MAIN[p_idx] << "\t" << SE_MAIN[p_idx] << "\t" << "NA" << "\t" << STAT_INT[p_idx] << "\t" << p_i << "\t" << "NA" << "\n";
            }
            p_idx++;
           }
         }
          
         npbidx = 0; 
         snp_skip.zeros();
         G.reshape(n, npb);
          
       }
      if((m+1) % 100000 == 0) {writefile << flush;}
     }
     if(end % 100000 != 0) {writefile << flush;}
     (void)ret;
     delete [] tmpout;
     delete [] snpID;
     delete [] rsID;
     delete [] chrStr;
     delete [] allele0;
     delete [] allele1;
     writefile.close();
     writefile.clear();
     fclose(fp);
     return wrap(compute_time);
    } 
     catch( std::exception &ex ) {
     forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;
  }  
  
  SEXP glmm_gei_bgen11(SEXP res_in, SEXP nullObj_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, 
                       SEXP ei_in, SEXP qi_in, SEXP isNullP_in, SEXP isNullEC_in,
                       SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
    try{
      
      bool isNullP = Rcpp::as<bool>(isNullP_in);
      bool isNullEC = Rcpp::as<bool>(isNullEC_in);
      int ei = Rcpp::as<int>(ei_in);
      int qi = Rcpp::as<int>(qi_in);
      
      arma::mat P;
      arma::sp_mat Sigma_i;
      arma::sp_mat Sigma_iX;
      arma::sp_mat cov;
      
      Rcpp::List null_obj(nullObj_in);
      
      arma::mat E = as<arma::mat>(null_obj["E"]);
      arma::mat EC;
      if (!isNullEC) {
        EC = as<arma::mat>(null_obj["EC"]);
      } else {
        (void)EC;
      }
      
      if (!isNullP) {
        arma::mat P = as<arma::mat>(null_obj["P"]);
        (void)Sigma_i; (void)Sigma_iX; (void)cov;
        
      } else {
        Sigma_i = as<arma::sp_mat>(null_obj["Sigma_i"]);
        Sigma_iX = as<arma::sp_mat>(null_obj["Sigma_iX"]);
        cov = as<arma::sp_mat>(null_obj["cov"]);
        (void)P;
      }
      
      
      Rcpp::NumericVector res_r(res_in);
      arma::vec res(res_r.begin(), res_r.size(), false);
      const char center = Rcpp::as<char>(center_in);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const char miss_method = Rcpp::as<char>(miss_method_in);
      
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      string outfile  = Rcpp::as<string>(outfile_in);
      
      const size_t npb = Rcpp::as<size_t>(nperbatch_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n), snp_skip = zeros<uvec>(npb);
      mat G(n, npb);
      string* tmpout = new string[npb];
      vector <string> biminfo;
      double gmean, geno, gmax, gmin;
      const double tol = 1e-5;
      size_t ncount, nmiss, npbidx = 0;
      double compute_time = 0.0;
      ofstream writefile(outfile.c_str(), ofstream::out);
      
      
      
      uint maxLA = 65536;
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];
      uint Nbgen = Rcpp::as<uint>(nbgen_in);
      uint compression = Rcpp::as<uint>(compression_in);
      std::vector<uchar> zBuf11;
      std::vector<uchar> shortBuf11;
      uLongf destLen1 = 6 * Nbgen;
      zBuf11.resize(destLen1);
      shortBuf11.resize(destLen1);
      
      uint begin = Rcpp::as<uint>(begin_in);
      uint end   = Rcpp::as<uint>(end_in);
      long long unsigned int byte = Rcpp::as<long long unsigned int>(pos_in);
      FILE* fp = fopen(bgenfile.c_str(), "rb");
      fseek(fp, byte, SEEK_SET);
      
      bool isMultiThread = Rcpp::as<bool>(isMultiThread_in);
      
      if (!isMultiThread){
        writefile << "SNP\tRSID\tCHR\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
      }
      
      int ret;
      
      for (uint m = begin; m < end; m++) {
        stringstream writeout;
        
        uint Nprob; ret = fread(&Nprob, 4, 1, fp); 
        if (Nprob != Nbgen) {
          Rcout << "Error reading bgen file: Number of samples with genotype probabilities does not match number of samples in BGEN file. \n"; return R_NilValue;
        }
        ushort LS; ret = fread(&LS, 2, 1, fp);
        ret = fread(snpID, 1, LS, fp); snpID[LS] = '\0';
        
        ushort LR; ret = fread(&LR, 2, 1, fp);
        ret = fread(rsID, 1, LR, fp);  rsID[LR] = '\0';
        
        ushort LC; ret = fread(&LC, 2, 1, fp);
        ret = fread(chrStr, 1, LC, fp); chrStr[LC] = '\0';
        
        uint physpos; ret = fread(&physpos, 4, 1, fp);
        string physpos_tmp = to_string(physpos);
        
        uint32_t LA; 
        ret = fread(&LA, 4, 1, fp);
        ret = fread(allele1, 1, LA, fp);
        allele1[LA] = '\0';
        
        uint32_t LB; 
        ret = fread(&LB, 4, 1, fp);
        ret = fread(allele0, 1, LB, fp); 
        allele0[LB] = '\0';
        
        
        uint16_t* probs_start;
        if (compression == 1) {
          uint cLen; 
          ret = fread(&cLen, 4, 1, fp);
          ret = fread(&zBuf11[0], 1, cLen, fp);
                        
          uLongf destLen = 6*Nbgen;
          if (uncompress(&shortBuf11[0], &destLen, &zBuf11[0], cLen) != Z_OK || destLen != (destLen1)) {
            Rcout << "Error reading bgen file: Decompressing variant block failed with zLib. \n"; return R_NilValue;
          }
          
          probs_start = reinterpret_cast<uint16_t*>(&shortBuf11[0]);
        }
        else {
          ret = fread(&zBuf11[0], 1, destLen1, fp);
          probs_start = reinterpret_cast<uint16_t*>(&zBuf11[0]);
        }
        
        const double scale = 1.0 / 32768;
        
        
        gmean=0.0;
        gmax=-100.0;
        gmin=100.0;
        nmiss=0;
        ncount = 0;
        for (size_t i = 0; i < Nbgen; i++) {
          if (select[ncount] > 0) {
            
            double p11 = probs_start[3 * i] * scale;
            double p10 = probs_start[3 * i + 1] * scale;
            double p00 = probs_start[3 * i + 2] * scale;
            if (p11 == 0.0 && p10 == 0.0 && p00 == 0.0){
              gmiss[select[ncount]-1] = 1;
              nmiss++;
              
            } else {
              double pTot = p11 + p10 + p00;
              geno = (2 * p00 + p10) / pTot;
              gmiss[select[ncount]-1] = 0;
              g[select[ncount]-1] = geno;
              gmean += geno;
              if (geno > gmax) { gmax = geno; }
              if (geno < gmin) { gmin = geno; }
            } 
            
          }
          ncount++;
        }
        
        gmean/=(double)(n-nmiss);
        
        for (size_t j=0; j<n; ++j) {
          if (gmiss[j]==1) {
            g[j] = gmean;
            if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
          }
          if (center=='c') {
            g[j] -= gmean;
          }
        }
        
        gmean /= 2.0; // convert mean to allele freq
        writeout << snpID << "\t" << rsID << "\t" << chrStr << "\t" << physpos_tmp << "\t" << allele1 << "\t" << allele0 << "\t" << (n-nmiss) << "\t" << (nmiss/n) << "\t";
        tmpout[npbidx] = writeout.str();
        writeout.clear();
        if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
          snp_skip[npbidx] = 1;
        } else {
          G.col(npbidx) = g;
        }
        npbidx++;
        
        
        vec V_i;
        mat IV_U;
        mat IV_V_i;
        vec BETA_MAIN;
        vec STAT_MAIN;
        vec STAT_INT;
        vec SE_MAIN;
        if((m+1 == end) || (npbidx == npb)) {
          
          
          uvec snp_idx = find(snp_skip == 0);
          if (npbidx != npb) {
            G.reshape(n, npbidx);
            snp_idx = snp_idx.rows(0,npbidx-1);
          }
          G = G.cols(snp_idx);
          
          if (G.n_cols != 0) {
            int ng = G.n_cols;
            mat K;
            for (int ei_idx = 0; ei_idx < ei; ei_idx++) {
              mat tmpK = G.each_col() % E.col(ei_idx);
              K = join_horiz(K, tmpK);
            } 
            if (!isNullEC) {
              uvec g_idx = linspace<uvec>(0, ng-1, ng);
              for (int qi_idx = 0; qi_idx < qi; qi_idx++) {
                mat tmpG = G.cols(g_idx);
                tmpG = tmpG.each_col() % EC.col(qi_idx);
                G = join_horiz(G, tmpG);
              }
            }
            
            vec U = G.t() * res;
            sp_mat Gsp(G);
            sp_mat PG;
            if (!isNullP) {
              PG = P.t() * G;
            } else {
              sp_mat GSigma_iX =  Gsp.t() * Sigma_iX;
              PG = (Sigma_i.t() * Gsp) - (Sigma_iX * (GSigma_iX * cov.t()).t());
            }
            mat GPG = (G.t() * PG)  % kron(ones(1+qi, 1+qi), mat(ng, ng, fill::eye));
            mat GPG_i = inv(GPG);
            
            if (isNullEC) {
              V_i = diagvec(GPG_i);
            } else {
              V_i = diagvec(GPG);
              V_i = V_i.rows(0, ng-1);
              double eps = DBL_EPSILON;
              for (size_t tmp_v = 0; tmp_v < V_i.size(); tmp_v++) {
                if (V_i[tmp_v] < eps) {
                  V_i[tmp_v] = 0;
                } else {
                  V_i[tmp_v] =  1 / V_i[tmp_v];
                }
              }
            }
            
            BETA_MAIN = V_i % U.rows(0,ng-1);
            SE_MAIN = sqrt(V_i);
            STAT_MAIN = BETA_MAIN % U.rows(0,ng-1);
            mat KPK;
            
            if (!isNullP) {
              KPK = K.t() * (P.t(), K);
            } else {
              mat KSigma_iX = K.t() * Sigma_iX;
              KPK = (K.t() * (Sigma_i.t() * K)) - (KSigma_iX * (KSigma_iX * cov.t()).t());
            }
            
            KPK = KPK % kron(ones(ei, ei), mat(ng, ng, fill::eye));
            mat KPG = (K.t() * PG) % kron(ones(ei, 1+qi), mat(ng, ng, fill::eye));
            if (isNullEC) {
              mat IV_U_kron = kron(ones(ei,1), mat(ng, ng, fill::eye));
              IV_U = IV_U_kron.each_col() % ((K.t() * res) - (KPG * BETA_MAIN));
              IV_V_i = inv(KPK - (KPG * (KPG * GPG_i.t()).t()));
              STAT_INT = diagvec(IV_U.t() * (IV_V_i.t() * IV_U));
              
            } else {
              mat IV_U_kron = kron(ones(ei,1), mat(ng, ng, fill::eye));
              IV_U = IV_U_kron.each_col() % ((K.t() * res) - (KPG * (GPG_i.t() * U)));
              IV_V_i = inv(KPK - (KPG * (KPG * GPG_i.t()).t()));
              STAT_INT = diagvec(IV_U.t() * (IV_V_i.t() * IV_U));
            }
            
          }
          
          int p_idx = 0;
          for(size_t j=0; j<npbidx; ++j) {
            if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
              writefile << tmpout[j] << "NA\tNA\tNA\tNA\tNA\tNA\n";
            } 
            else {
              if (V_i[j] > 0) {
                double p_m = Rf_pchisq(STAT_MAIN[p_idx], 1, 0, 0);
                double p_i = Rf_pchisq(STAT_INT[p_idx], ei, 0, 0);
                double p_j = Rf_pchisq(STAT_MAIN[p_idx] + STAT_INT[p_idx], 1+ei, 0, 0);
                writefile << tmpout[j] << BETA_MAIN[p_idx] << "\t" << SE_MAIN[p_idx] << "\t" << p_m << "\t" << STAT_INT[p_idx] << "\t" << p_i << "\t" << p_j << "\n";
              } 
              else {
                double p_i = Rf_pchisq(STAT_INT[j], ei, 0, 0);
                writefile << tmpout[j] << BETA_MAIN[p_idx] << "\t" << SE_MAIN[p_idx] << "\t" << "NA" << "\t" << STAT_INT[p_idx] << "\t" << p_i << "\t" << "NA" << "\n";
              }
              p_idx++;
            }
          }
          
          npbidx = 0; 
          snp_skip.zeros();
          G.reshape(n, npb);
          
        }
        if((m+1) % 100000 == 0) {writefile << flush;}
      }
      
      
      if(end % 100000 != 0) {writefile << flush;}
      (void)ret;
      delete [] tmpout;
      delete [] snpID;
      delete [] rsID;
      delete [] chrStr;
      delete [] allele0;
      delete [] allele1;
      writefile.close();
      writefile.clear();
      fclose(fp);
      return wrap(compute_time);
    } 
    catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;
  }  
} // end of extern "C"
