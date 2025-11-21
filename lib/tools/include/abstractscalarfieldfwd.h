#ifndef ABSTRACTSCALARFIELDFWD_H
#define ABSTRACTSCALARFIELDFWD_H

#include <memory>

namespace IPP
{

class AbstractScalarField;
typedef std::shared_ptr<AbstractScalarField> AbstractScalarFieldPtr;

class ConstAbstractScalarField;
typedef std::shared_ptr<ConstAbstractScalarField> ConstAbstractScalarFieldPtr;

}

#endif // ABSTRACTSCALARFIELDFWD_H
