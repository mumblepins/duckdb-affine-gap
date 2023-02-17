#define DUCKDB_EXTENSION_MAIN
#include "affine_gap_extension.hpp"
#include "duckdb/function/scalar/string_functions.hpp"
#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector_operations/septenary_executor.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"

#include "duckdb/function/scalar_function.hpp"

#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

namespace duckdb {

static double AffineGapDistance(const string_t &s1, const string_t &s2,
                                double matchWeight = 1,
                                double mismatchWeight = 11,
                                double gapWeight = 10, double spaceWeight = 7,
                                double abbreviationScale = .125) {

  auto s1_len = s1.GetSize();
  auto s2_len = s2.GetSize();
  auto s1_str = s1.GetDataUnsafe();
  auto s2_str = s2.GetDataUnsafe();
  if (s1_len == s2_len && s1_str == s2_str && matchWeight <= mismatchWeight &&
      matchWeight <= gapWeight) {
    return matchWeight * s1_len;
  }

  if (s1_len < s2_len) {
    std::swap(s1_len, s2_len);
    std::swap(s1_str, s2_str);
  }
  std::vector<double> D(s1_len + 1, 0);
  std::vector<double> V_current(s1_len + 1, 0);
  std::vector<double> V_previous(s1_len + 1, 0);

  double distance;
  char char1;
  char char2;
  double I, M;
  V_current[0] = 0;
  for (idx_t j = 1; j <= s1_len; j++) {
    V_current[j] = gapWeight + spaceWeight * j;
    D[j] = UINT64_MAX;
  }

  for (idx_t i = 1; i <= s2_len; i++) {
    char2 = s2_str[i - 1];
    for (idx_t j = 0; j <= s1_len; j++) {
      V_previous[j] = V_current[j];
    }

    // # Base conditions
    // # V(i,0) = gapWeight + spaceWeight * i
    // # I(i,0) = Infinity
    V_current[0] = gapWeight + spaceWeight * i;
    I = INT64_MAX;
    for (idx_t j = 1; j <= s1_len; j++) {
      char1 = s1_str[j - 1];
      // I(i,j) is the edit distance if the jth character of string 1
      // was inserted into string 2.
      //
      // I(i,j) = min(I(i,j-1), V(i,j-1) + gapWeight) + spaceWeight
      if (j <= s2_len) {
        I = std::min(I, V_current[j - 1] + gapWeight) + spaceWeight;
      } else {
        // Pay less for abbreviations
        // i.e. 'spago (los angeles) to 'spago'
        I = (std::min(I, V_current[j - 1] + gapWeight * abbreviationScale) +
             spaceWeight * abbreviationScale);
      }

      // D(i,j) is the edit distance if the ith character of string 2
      // was deleted from string 1
      //
      // D(i,j) = min((i-1,j), V(i-1,j) + gapWeight) + spaceWeight
      D[j] = std::min(D[j], V_previous[j] + gapWeight) + spaceWeight;

      // M(i,j) is the edit distance if the ith and jth characters
      // match or mismatch
      //
      // M(i,j) = V(i-1,j-1) + (matchWeight | misMatchWeight)
      if (char2 == char1) {
        M = V_previous[j - 1] + matchWeight;
      } else {
        M = V_previous[j - 1] + mismatchWeight;
      }

      // V(i,j) is the minimum edit distance
      //
      // V(i,j) = min(E(i,j), F(i,j), G(i,j))
      V_current[j] = std::min(std::min(I, D[j]), M);
    }
  }

  return V_current[s1_len];
}

static double
NormalizedAffineGapDistance(const string_t &s1, const string_t &s2,
                            double matchWeight = 1, double mismatchWeight = 11,
                            double gapWeight = 10, double spaceWeight = 7,
                            double abbreviationScale = .125) {
  idx_t len1 = s1.GetSize();
  idx_t len2 = s2.GetSize();
  if (len1 == 0 && len2 == 0) {
    return 0.5;
  }
  double distance =
      AffineGapDistance(s1, s2, matchWeight, mismatchWeight, gapWeight,
                        spaceWeight, abbreviationScale);
  return distance / (len1 + len2);
}

static double AffineGapDistanceScalarFunction(
    Vector &result, const string_t s1, const string_t s2,
    double matchWeight = 1, double mismatchWeight = 11, double gapWeight = 10,
    double spaceWeight = 7, double abbreviationScale = .125) {
  return (double)AffineGapDistance(s1, s2, matchWeight, mismatchWeight,
                                   gapWeight, spaceWeight, abbreviationScale);
}

static double NormalizedAffineGapDistanceScalarFunction(
    Vector &result, const string_t s1, const string_t s2,
    double matchWeight = 1, double mismatchWeight = 11, double gapWeight = 10,
    double spaceWeight = 7, double abbreviationScale = .125) {
  return (double)NormalizedAffineGapDistance(s1, s2, matchWeight,
                                             mismatchWeight, gapWeight,
                                             spaceWeight, abbreviationScale);
}

static void AffineGapFunction(DataChunk &args, ExpressionState &state,
                              Vector &result) {
  D_ASSERT(args.ColumnCount() == 2 || args.ColumnCount() == 7);
  auto &s1_vec = args.data[0];
  auto &s2_vec = args.data[1];
  if (args.ColumnCount() == 2) {
    BinaryExecutor::Execute<string_t, string_t, double>(
        s1_vec, s2_vec, result, args.size(),

        [&](string_t s1, string_t s2) {
          return AffineGapDistanceScalarFunction(result, s1, s2);
        });
  } else {
    SeptenaryExecutor::Execute<string_t, string_t, double, double, double,
                               double, double, double>(
        args, result,
        [&](string_t s1, string_t s2, double mw, double mmw, double gw,
            double sw, double as) {
          return AffineGapDistanceScalarFunction(result, s1, s2, mw, mmw, gw,
                                                 sw, as);
        });
  }
}

static void NormalizedAffineGapFunction(DataChunk &args, ExpressionState &state,
                                        Vector &result) {
  D_ASSERT(args.ColumnCount() == 2 || args.ColumnCount() == 7);
  auto &s1_vec = args.data[0];
  auto &s2_vec = args.data[1];
  if (args.ColumnCount() == 2) {
    BinaryExecutor::Execute<string_t, string_t, double>(
        s1_vec, s2_vec, result, args.size(),

        [&](string_t s1, string_t s2) {
          return NormalizedAffineGapDistanceScalarFunction(result, s1, s2);
        });
  } else {
    SeptenaryExecutor::Execute<string_t, string_t, double, double, double,
                               double, double, double>(
        args, result,
        [&](string_t s1, string_t s2, double mw, double mmw, double gw,
            double sw, double as) {
          return NormalizedAffineGapDistanceScalarFunction(result, s1, s2, mw,
                                                           mmw, gw, sw, as);
        });
  }
}

static void LoadInternal(DatabaseInstance &instance) {
  Connection con(instance);
  con.BeginTransaction();

  auto &catalog = Catalog::GetSystemCatalog(*con.context);

  ScalarFunctionSet affine_gap("affine_gap");
  affine_gap.AddFunction(
      ScalarFunction("affine_gap", {LogicalType::VARCHAR, LogicalType::VARCHAR},
                     LogicalType::DOUBLE, AffineGapFunction));
  affine_gap.AddFunction(ScalarFunction(
      "affine_gap",
      {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE,
       LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE,
       LogicalType::DOUBLE},
      LogicalType::DOUBLE, AffineGapFunction));
  CreateScalarFunctionInfo affine_gap_fun_info(affine_gap);

  affine_gap_fun_info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
  catalog.CreateFunction(*con.context, &affine_gap_fun_info);

  ScalarFunctionSet normalized_affine_gap("normalized_affine_gap");
  normalized_affine_gap.AddFunction(ScalarFunction(
      "normalized_affine_gap", {LogicalType::VARCHAR, LogicalType::VARCHAR},
      LogicalType::DOUBLE, NormalizedAffineGapFunction));
  normalized_affine_gap.AddFunction(ScalarFunction(
      "normalized_affine_gap",
      {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE,
       LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE,
       LogicalType::DOUBLE},
      LogicalType::DOUBLE, NormalizedAffineGapFunction));
  CreateScalarFunctionInfo normalized_affine_gap_fun_info(
      normalized_affine_gap);

  normalized_affine_gap_fun_info.on_conflict =
      OnCreateConflict::ALTER_ON_CONFLICT;
  catalog.CreateFunction(*con.context, &normalized_affine_gap_fun_info);
  con.Commit();
}

void Affine_gapExtension::Load(DuckDB &db) { LoadInternal(*db.instance); }
std::string Affine_gapExtension::Name() { return "affine_gap"; }

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void affine_gap_init(duckdb::DatabaseInstance &db) {
  LoadInternal(db);
}

DUCKDB_EXTENSION_API const char *affine_gap_version() {
  return duckdb::DuckDB::LibraryVersion();
}
}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
