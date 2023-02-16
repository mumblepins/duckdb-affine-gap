#pragma once

#include "duckdb.hpp"

namespace duckdb {

class Affine_gapExtension : public Extension {
public:
  void Load(DuckDB &db) override;
  std::string Name() override;
};

} // namespace duckdb
