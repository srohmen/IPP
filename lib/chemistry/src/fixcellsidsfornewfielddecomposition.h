#ifndef FIXCELLSIDSFORNEWFIELDDECOMPOSITION_H
#define FIXCELLSIDSFORNEWFIELDDECOMPOSITION_H

#include <string>

namespace IPP
{

class FieldDecomposition;

namespace FixCellsIDsForNewFieldDecomposition
{

void fix(const FieldDecomposition& oldDecomp,
         const FieldDecomposition& newDecomp,
         std::istream& input,
         std::ostream& output);

}

}

#endif // FIXCELLSIDSFORNEWFIELDDECOMPOSITION_H
