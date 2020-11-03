#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP glmm_gei_bgen11(SEXP res_in, SEXP P_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, 
                     SEXP ei_in, SEXP qi_in, SEXP isNullP_in, SEXP isNullEC_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

SEXP glmm_gei_bgen13(SEXP res_in, SEXP P_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, 
                     SEXP ei_in, SEXP qi_in, SEXP isNullP_in, SEXP isNullEC_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

SEXP bgenHeader(SEXP bgenfile_in);

SEXP getVariantPos(SEXP bgenfile_in, SEXP offset_in, SEXP mbgen_in, SEXP nbgen_in, SEXP compression_in, SEXP layout_in, SEXP cores_in);

SEXP bgenVariantInfo(SEXP bgenfile_in, SEXP offset_in, SEXP mbgen_in, SEXP nbgen_in, SEXP layout_in, SEXP compression_in);

SEXP magee_bgen13(SEXP res_in, SEXP nullObj_in, SEXP bgenfile_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, 
                  SEXP isNullP_in, SEXP n_groups_in, SEXP mf_in, SEXP mv_in, SEXP if_in, SEXP iv_in, SEXP jd_in, SEXP jf_in, SEXP jv_in,
                  SEXP groupInfo_in, SEXP groupStart_in, SEXP groupEnd_in, SEXP weight_in, SEXP method_in, SEXP beta1_in, SEXP beta2_in,
                  SEXP select_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in);

static const R_CallMethodDef R_CallDef[]  = {
  {"magee_bgen13", (DL_FUNC) &magee_bgen13, 27},
  {"glmm_gei_bgen11", (DL_FUNC) &glmm_gei_bgen11, 23},
  {"glmm_gei_bgen13", (DL_FUNC) &glmm_gei_bgen13, 23},
  {"bgenHeader", (DL_FUNC) &bgenHeader, 1},
  {"getVariantPos", (DL_FUNC) &getVariantPos, 7},
  {"bgenVariantInfo", (DL_FUNC) &bgenVariantInfo, 6},
  {NULL, NULL, 0}
};

void R_init_GMMAT(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
