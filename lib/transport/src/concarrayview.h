#ifndef CONCARRAYVIEW_H
#define CONCARRAYVIEW_H

#include <cstddef>
#include "array_view.h"
#include "simulationexchangedata.h"
#include "ippbox.h"
#include "fielddecomposition.h"

namespace IPP
{

namespace ConcArrayViewTool
{

template<size_t dim, typename Container>
av::array_view<typename Container::value_type, dim+1> makeConcView(const size_t nComps,
                                                                   const size_t nx,
                                                                   const size_t ny,
                                                                   const size_t nz,
                                                                   Container& container)
{
    const av::bounds<dim+1> bounds = av::MakeBounds<dim>::get((ptrdiff_t)nComps,
                                                              (ptrdiff_t)nx,
                                                              (ptrdiff_t)ny,
                                                              (ptrdiff_t)nz);
    using T = typename Container::value_type;
    av::array_view<T, dim+1> view(container, bounds);
    return view;
}

template<size_t dim, typename Size, typename Container>
av::array_view<typename Container::value_type, dim+1> makeConcView(const size_t nComps,
                                                                   const Size& size,
                                                                   Container& container)
{
    return makeConcView<dim>(nComps, size[0], size[1], size[2], container);
}



template<size_t dim, typename Container>
av::array_view<const typename Container::value_type, dim+1> makeConcView(const size_t nComps,
                                                                         const size_t nx,
                                                                         const size_t ny,
                                                                         const size_t nz,
                                                                         const Container& container)
{
    const av::bounds<dim+1> bounds = av::MakeBounds<dim>::get((ptrdiff_t)nComps,
                                                              (ptrdiff_t)nx,
                                                              (ptrdiff_t)ny,
                                                              (ptrdiff_t)nz);
    av::array_view<const typename Container::value_type, dim+1> view(container, bounds);
    return view;
}

template<size_t dim, typename Size, typename Container>
av::array_view<const typename Container::value_type, dim+1> makeConcView(const size_t nComps,
                                                                         const Size& size,
                                                                         const Container& container)
{
    return makeConcView<dim>(nComps, size[0], size[1], size[2], container);
}

template<size_t dim, typename Container>
av::array_view<typename Container::value_type, dim+1> makeConcView(const SimulationExchangeData& data,
                                                                   Container& concRaw)
{
    const FieldDecomposition* decomp = data.getDecomp();
    const IPPBox3DLong& ownField = decomp->getOwnDomain();
    const IPPVector3DLong ownSize = ownField.getSize();
    const size_t nComps = data.getCompNames().size();

    const av::array_view<typename Container::value_type, dim+1> concView = makeConcView<dim>(nComps, ownSize, concRaw);

    return concView;
}

template<size_t dim, typename Container>
av::array_view<const typename Container::value_type, dim+1> makeConcView(const SimulationExchangeData& data,
                                                                         const Container& concRaw)
{
    const FieldDecomposition* decomp = data.getDecomp();
    const IPPBox3DLong& ownField = decomp->getOwnDomain();
    const IPPVector3DLong ownSize = ownField.getSize();
    const size_t nComps = data.getCompNames().size();

    const av::array_view<const typename Container::value_type, dim+1> concView = makeConcView<dim>(nComps, ownSize, concRaw);

    return concView;
}

} // end of ConcArrayViewTool

} // end of IPP

#endif // CONCARRAYVIEW_H
