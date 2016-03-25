#include "def.hpp"

struct DOFobject
{
  DOFobject():
    dof_id(DOFidNotSet)
  {}

  DOFid_T dof_id;
};
