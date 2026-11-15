#pragma once

#include "implicits.h"
#include "sdf_node.h"

#include <stdexcept>
#include <utility>

// Adapts a construction tree to TinyMesh's existing polygonizer and normals.
class SdfField final : public AnalyticScalarField
{
public:
  explicit SdfField(sdf::NodePtr root) : root_(std::move(root))
  {
    if (!root_)
      throw std::invalid_argument("An SDF field requires a root node");
  }

  double Value(const Vector& point) const override
  {
    return root_->Value(point);
  }

private:
  sdf::NodePtr root_;
};
