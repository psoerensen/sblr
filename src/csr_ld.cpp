#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>

struct CSR {
 std::vector<int> indptr;
 std::vector<int> indices;
 std::vector<float> values;
};

CSR readLD_to_CSR_cpp(
  const std::string& filename,
  int mchr,
  int msize,
  double ld_threshold,   // renamed for clarity
  bool onebased = false
) {
 std::ifstream file(filename, std::ios::binary);
 if (!file.is_open()) {
  throw std::runtime_error("Cannot open LD file.");
 }

 const int window = 2 * msize + 1;

 // ----------------------------
 // Validate file size (IMPORTANT)
 // ----------------------------
 file.seekg(0, std::ios::end);
 size_t filesize = file.tellg();
 file.seekg(0, std::ios::beg);

 size_t expected = static_cast<size_t>(mchr) * window * sizeof(float);

 if (filesize != expected) {
  throw std::runtime_error("LD file size mismatch.");
 }

 // ----------------------------
 // Allocate CSR
 // ----------------------------
 CSR csr;
 csr.indptr.resize(mchr + 1);

 // Reserve enough memory to avoid reallocations
 csr.indices.reserve(static_cast<size_t>(mchr) * window);
 csr.values.reserve(static_cast<size_t>(mchr) * window);

 size_t nnz = 0;
 csr.indptr[0] = 0;

 std::vector<float> buffer(window);

 // ----------------------------
 // Read file column-wise
 // ----------------------------
 for (int i = 0; i < mchr; ++i) {

  file.read(reinterpret_cast<char*>(buffer.data()),
            window * sizeof(float));

  if (!file) {
   throw std::runtime_error("Error reading LD file.");
  }

  for (int k = 0; k < window; ++k) {

   float ld = buffer[k];

   // force diagonal = 1
   if (k == msize) ld = 1.0f;

   // skip invalid values
   if (!std::isfinite(ld)) continue;

   // threshold on absolute LD (better than r^2)
   if (std::abs(ld) > ld_threshold) {

    int j = i + k - msize;

    if (j < 0 || j >= mchr) continue;

    csr.indices.push_back(onebased ? j + 1 : j);
    csr.values.push_back(ld);
    ++nnz;
   }
  }

  csr.indptr[i + 1] = static_cast<int>(nnz);
 }

 file.close();
 return csr;
}

extern "C" {
#include <R.h>
#include <Rinternals.h>
}

// [[Rcpp::export]]
SEXP readLD_to_CSR_R(
  SEXP filenameSEXP,
  SEXP mchrSEXP,
  SEXP msizeSEXP,
  SEXP thresholdSEXP,
  SEXP onebasedSEXP
) {
 const char* filename = CHAR(STRING_ELT(filenameSEXP, 0));
 int mchr = INTEGER(mchrSEXP)[0];
 int msize = INTEGER(msizeSEXP)[0];
 double threshold = REAL(thresholdSEXP)[0];
 bool onebased = LOGICAL(onebasedSEXP)[0];

 CSR csr = readLD_to_CSR_cpp(filename, mchr, msize, threshold, onebased);

 // ----------------------------
 // Convert to R objects
 // ----------------------------
 SEXP indptr = PROTECT(allocVector(INTSXP, csr.indptr.size()));
 SEXP indices = PROTECT(allocVector(INTSXP, csr.indices.size()));
 SEXP values  = PROTECT(allocVector(REALSXP, csr.values.size()));

 for (size_t i = 0; i < csr.indptr.size(); ++i)
  INTEGER(indptr)[i] = csr.indptr[i];

 for (size_t i = 0; i < csr.indices.size(); ++i)
  INTEGER(indices)[i] = csr.indices[i];

 for (size_t i = 0; i < csr.values.size(); ++i)
  REAL(values)[i] = csr.values[i];

 SEXP result = PROTECT(allocVector(VECSXP, 3));
 SET_VECTOR_ELT(result, 0, indptr);
 SET_VECTOR_ELT(result, 1, indices);
 SET_VECTOR_ELT(result, 2, values);

 SEXP names = PROTECT(allocVector(STRSXP, 3));
 SET_STRING_ELT(names, 0, mkChar("indptr"));
 SET_STRING_ELT(names, 1, mkChar("indices"));
 SET_STRING_ELT(names, 2, mkChar("values"));

 setAttrib(result, R_NamesSymbol, names);

 UNPROTECT(5);
 return result;
}
