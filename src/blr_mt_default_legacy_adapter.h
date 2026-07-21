#ifndef SBLR_BLR_MT_DEFAULT_LEGACY_ADAPTER_H
#define SBLR_BLR_MT_DEFAULT_LEGACY_ADAPTER_H

#include "blr_mt_default_types.h"

#include <vector>

namespace sblr {
namespace mt {

using MtDefaultLegacyResult=
 std::vector<std::vector<std::vector<double>>>;

inline MtDefaultLegacyResult make_mt_default_legacy_result(
 const MtDefaultFinalResult& value,
 const std::vector<std::vector<double>>& wy,
 int nit,
 int nburn
) {
 const int nt=value.nt;
 const int m=value.m;
 const int nmodels=value.nmodels;
 MtDefaultLegacyResult result(20);
 for (auto& field : result) field.resize(nt);
 for (int t=0; t<nt; ++t) {
  for (int field=0; field<=6; ++field) result[field][t].resize(m);
  for (int field=7; field<=9; ++field) result[field][t].resize(nit+nburn);
  for (int field=10; field<=15; ++field) result[field][t].resize(nt);
  result[16][t].resize(nmodels);
  result[17][t].resize(nmodels);
  result[18][t].resize(4);
  result[19][t].resize(2);
 }
 for (int t=0; t<nt; ++t) {
  for (int i=0; i<m; ++i) {
   result[0][t][i]=value.bm[t][i];
   result[1][t][i]=value.dm[t][i];
   result[2][t][i]=wy[t][i];
   result[3][t][i]=value.r[t][i];
   result[4][t][i]=value.b[t][i];
   result[5][t][i]=value.d[t][i];
   result[6][t][i]=value.marker_order[i];
  }
  for (int i=0; i<nit+nburn; ++i) {
   result[7][t][i]=value.vbs[t][i];
   result[8][t][i]=value.vgs[t][i];
   result[9][t][i]=value.ves[t][i];
  }
  for (int t2=0; t2<nt; ++t2) {
   result[10][t][t2]=value.covb[t][t2];
   result[11][t][t2]=value.covg[t][t2];
   result[12][t][t2]=value.cove[t][t2];
   result[13][t][t2]=value.vb(t,t2);
   result[14][t][t2]=value.vg(t,t2);
   result[15][t][t2]=value.ve(t,t2);
  }
  for (int i=0; i<nmodels; ++i) {
   result[16][t][i]=value.pi_final[i];
   result[17][t][i]=value.pi_mean[i];
  }
  for (int i=0; i<4; ++i) result[18][t][i]=value.pitrait[t][i];
  for (int i=0; i<2; ++i) result[19][t][i]=value.pimarker[i];
 }
 return result;
}

}  // namespace mt
}  // namespace sblr

#endif
