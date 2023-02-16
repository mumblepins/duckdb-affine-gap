#define DUCKDB_EXTENSION_MAIN

#include "affine_gap_extension.hpp"
#include "duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/scalar_function.hpp"

#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

namespace duckdb {

inline void Affine_gapScalarFun(DataChunk &args, ExpressionState &state,
                                Vector &result) {
  auto &name_vector = args.data[0];
  UnaryExecutor::Execute<string_t, string_t>(
      name_vector, result, args.size(), [&](string_t name) {
        return StringVector::AddString(result,
                                       "Affine_gap " + name.GetString() + " 🐥");
        ;
      });
}

static void LoadInternal(DatabaseInstance &instance) {
  Connection con(instance);
  con.BeginTransaction();

  auto &catalog = Catalog::GetSystemCatalog(*con.context);

  CreateScalarFunctionInfo affine_gap_fun_info(
      ScalarFunction("affine_gap", {LogicalType::VARCHAR}, LogicalType::VARCHAR,
                     Affine_gapScalarFun));
  affine_gap_fun_info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
  catalog.CreateFunction(*con.context, &affine_gap_fun_info);
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
