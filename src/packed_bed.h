#ifndef PACKED_BED_H
#define PACKED_BED_H

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

static inline std::size_t round_up(std::size_t x, std::size_t a) {
 return ((x + a - 1) / a) * a;
}

static inline void* aligned_malloc64(std::size_t nbytes) {
 const std::size_t nbytes_aligned =
  round_up(nbytes, static_cast<std::size_t>(64));

#if defined(_MSC_VER)
 return _aligned_malloc(nbytes_aligned, 64);
#elif defined(_WIN32)
 return std::malloc(nbytes_aligned);
#else
 void* ptr = nullptr;
 if (posix_memalign(&ptr, 64, nbytes_aligned) != 0) return nullptr;
 return ptr;
#endif
}

static inline void aligned_free64(void* ptr) {
#if defined(_MSC_VER)
 _aligned_free(ptr);
#else
 std::free(ptr);
#endif
}

static inline uint8_t get_bed_code(const uint8_t* src, int sample0) {
 const int byte_index = sample0 >> 2;
 const int shift = (sample0 & 3) << 1;
 return static_cast<uint8_t>((src[byte_index] >> shift) & 3u);
}

static inline void set_bed_code(uint8_t* dst, int sample0, uint8_t code) {
 const int byte_index = sample0 >> 2;
 const int shift = (sample0 & 3) << 1;
 dst[byte_index] &= static_cast<uint8_t>(~(3u << shift));
 dst[byte_index] |= static_cast<uint8_t>((code & 3u) << shift);
}


struct PackedBedMatrix {
 int n = 0;
 int m = 0;
 std::size_t nbytes = 0;
 std::size_t stride = 0;
 uint8_t* data = nullptr;

 PackedBedMatrix() = default;

 PackedBedMatrix(int n_, int m_) : n(n_), m(m_) {
  nbytes = static_cast<std::size_t>((n + 3) / 4);
  stride = round_up(nbytes, static_cast<std::size_t>(64));

  const std::size_t total = static_cast<std::size_t>(m) * stride;
  data = static_cast<uint8_t*>(aligned_malloc64(total));
  if (!data) throw std::runtime_error("Failed to allocate packed BED matrix.");

  std::memset(data, 0, total);
 }

 ~PackedBedMatrix() {
  if (data) aligned_free64(data);
 }

 PackedBedMatrix(const PackedBedMatrix&) = delete;
 PackedBedMatrix& operator=(const PackedBedMatrix&) = delete;

 PackedBedMatrix(PackedBedMatrix&& other) noexcept
  : n(other.n),
    m(other.m),
    nbytes(other.nbytes),
    stride(other.stride),
    data(other.data) {
  other.n = 0;
  other.m = 0;
  other.nbytes = 0;
  other.stride = 0;
  other.data = nullptr;
 }

 PackedBedMatrix& operator=(PackedBedMatrix&& other) noexcept {
  if (this != &other) {
   if (data) aligned_free64(data);

   n = other.n;
   m = other.m;
   nbytes = other.nbytes;
   stride = other.stride;
   data = other.data;

   other.n = 0;
   other.m = 0;
   other.nbytes = 0;
   other.stride = 0;
   other.data = nullptr;
  }
  return *this;
 }

 uint8_t* row(int j) {
  return data + static_cast<std::size_t>(j) * stride;
 }

 const uint8_t* row(int j) const {
  return data + static_cast<std::size_t>(j) * stride;
 }
};

static PackedBedMatrix read_bedfiles_to_packed_matrix(
  const std::vector<std::string>& bed_files,
  int n_bed,
  const int* rows0,
  int n_rows,
  const std::vector<std::vector<int>>& cls_by_file
) {
 const bool use_rows = rows0 != nullptr;
 const int n_used = use_rows ? n_rows : n_bed;
 const std::size_t nbytes_bed =
  static_cast<std::size_t>((n_bed + 3) / 4);

 int m_total = 0;
 for (const auto& x : cls_by_file) {
  m_total += static_cast<int>(x.size());
 }

 PackedBedMatrix G(n_used, m_total);
 std::vector<uint8_t> src_row(use_rows ? nbytes_bed : 0);

 int out_i = 0;

 for (std::size_t f = 0; f < bed_files.size(); ++f) {
  FILE* fs = std::fopen(bed_files[f].c_str(), "rb");
  if (!fs) {
   throw std::runtime_error("Could not open BED file: " + bed_files[f]);
  }

  try {
   unsigned char header[3];

   if (std::fread(header, 1, 3, fs) != 3) {
    throw std::runtime_error("Could not read BED header: " + bed_files[f]);
   }

   if (header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
    throw std::runtime_error(
      "BED file is not SNP-major PLINK BED: " + bed_files[f]
    );
   }

   for (int cls : cls_by_file[f]) {
    const long long offset =
     3LL +
     static_cast<long long>(cls - 1) *
     static_cast<long long>(nbytes_bed);

    if (std::fseek(fs, offset, SEEK_SET) != 0) {
     throw std::runtime_error(
       "fseek failed while reading BED row: " + bed_files[f]
     );
    }

    uint8_t* dst = G.row(out_i);

    if (!use_rows) {
     const std::size_t got = std::fread(dst, 1, nbytes_bed, fs);
     if (got != nbytes_bed) {
      throw std::runtime_error(
        "fread failed while reading BED row: " + bed_files[f]
      );
     }
    } else {
     const std::size_t got = std::fread(src_row.data(), 1, nbytes_bed, fs);
     if (got != nbytes_bed) {
      throw std::runtime_error(
        "fread failed while reading BED row: " + bed_files[f]
      );
     }

     std::memset(dst, 0, G.nbytes);

     for (int j = 0; j < n_rows; ++j) {
      set_bed_code(dst, j, get_bed_code(src_row.data(), rows0[j]));
     }
    }

    ++out_i;
   }

   std::fclose(fs);
   fs = nullptr;
  } catch (...) {
   if (fs) std::fclose(fs);
   throw;
  }
 }

 return G;
}

static inline double code_to_dosage(uint8_t code) {
 switch (code & 3u) {
 case 0u: return 2.0;
 case 1u: return std::numeric_limits<double>::quiet_NaN();
 case 2u: return 1.0;
 default: return 0.0;
 }
}
static std::vector<double> packed_to_dense_vector(
  const PackedBedMatrix& G
) {
 std::vector<double> X(static_cast<std::size_t>(G.n) * G.m);

 for (int j = 0; j < G.m; ++j) {
  const uint8_t* row = G.row(j);

  for (int i = 0; i < G.n; ++i) {
   const uint8_t code = get_bed_code(row, i);
   X[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * G.n] =
    code_to_dosage(code);
  }
 }

 return X;
}


static std::vector<double> compute_af_from_packed(const PackedBedMatrix& G) {
 std::vector<double> af(static_cast<std::size_t>(G.m));

#pragma omp parallel for schedule(static)
 for (int j = 0; j < G.m; ++j) {
  const uint8_t* packed = G.row(j);
  double sum = 0.0;
  int n_obs = 0;

  for (int i = 0; i < G.n; ++i) {
   const uint8_t code = get_bed_code(packed, i);

   switch (code & 3u) {
   case 0u:
    sum += 2.0;
    ++n_obs;
    break;
   case 1u:
    break;
   case 2u:
    sum += 1.0;
    ++n_obs;
    break;
   default:
    ++n_obs;
   break;
   }
  }

  af[static_cast<std::size_t>(j)] =
   n_obs > 0 ? sum / (2.0 * static_cast<double>(n_obs)) : 0.0;
 }

 return af;
}

#endif
