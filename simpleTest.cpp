#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include "polyscope/polyscope.h"

using namespace DGtal;
using namespace Z3i;

int main()
{
  polyscope::init();


  polyscope::show();
  return 0;
}
