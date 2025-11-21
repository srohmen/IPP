#ifndef SETNTENSORTOSOURCE_H
#define SETNTENSORTOSOURCE_H

#include "plbtypededuction.h"

#include "domainiterator.h"
#include "array_view.h"
#include "permeabiltyflag.h"
#include "dotfindneighbors.h"

namespace IPP
{


template<typename T,
         typename U,
         template<typename V> class Descriptor,
         typename MaskType>
class SetNTensorToSource : public PlbTypeDeduction::
        GetBoxProcessingFunctionalXD<Descriptor<T>::d>::value
{
public:
    static constexpr size_t dim = Descriptor<T>::d;
    using Dot = typename PlbTypeDeduction::GetDotXD<dim>::value;

    SetNTensorToSource(const Dot& location,
                       av::array_view<U, dim+1>& postReacConc)
        : m_location(location)
        , m_postReacConc(postReacConc)
    {

    }

    virtual SetNTensorToSource<T, U, Descriptor, MaskType>* clone() const override
    {
        return new SetNTensorToSource<T, U, Descriptor, MaskType>(*this);
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const override
    {
        assert(modified.size() >= 2);
        modified[0] = plb::modif::nothing; // mask
        modified[1] = plb::modif::nothing; // source NTensor

        for(size_t i = 1; i < modified.size(); ++i)
        {
            modified[i] = plb::modif::staticVariables; // lattice per comp
        }
    }


private:

    using MaskField = typename PlbTypeDeduction::GetScalarField_XD<MaskType, dim>::value;
    using NTensorField = typename PlbTypeDeduction::GetNTensorField_XD<T, dim>::value;
    using Lattice = typename PlbTypeDeduction::GetBlockLattice_XD<T, Descriptor>::value;


    using BaseClass =
    typename PlbTypeDeduction::GetBoxProcessingFunctionalXD<dim>::value;

    using Traits = PlbTypeDeduction::function_traits<decltype(&BaseClass::processGenericBlocks)>;
    using Box = typename Traits::template argument<1>::type;
    using AtomicBlockVec = typename Traits::template argument<2>::type;


public:
    virtual void processGenericBlocks(Box domain, AtomicBlockVec fields) override
    {
        // 0 = permMask
        // 1 = source NTensor
        // 2 -> n-1 = lattices

        const constexpr size_t startLattices = 2;
        assert(fields.size() >= startLattices);

        const size_t nComps = fields.size() - startLattices;

        std::vector<Dot> offsetVec;
        offsetVec.reserve(nComps);

        const MaskField& permMask = *dynamic_cast<MaskField*>(fields[0]);
        const NTensorField& sourceTensorField = *dynamic_cast<NTensorField*>(fields[1]);

        const Dot maskLocation = permMask.getLocation();
        const Dot offsetSource = plb::computeRelativeDisplacement(permMask, sourceTensorField);

        for(size_t iComp = 0; iComp < nComps; ++iComp)
        {
            Lattice& lattice = *static_cast<Lattice*>(fields[iComp + startLattices]);
            const Dot offset = plb::computeRelativeDisplacement(permMask, lattice);
            offsetVec.push_back(offset);
        }


        const auto end = DataAccess::end(domain);
        for(auto it = DataAccess::begin(domain); it < end; ++it)
        {
            const Dot& pos = *it;
            const MaskType& maskVal = DataAccess::get(permMask, pos);

            if(maskVal & PF_isPermeable)
            {
                const Dot absPos = pos + maskLocation - m_location;
                const auto neighbors = DotFindNeighbors::find(pos);

                for(const auto& neighPos : neighbors.dots)
                {
                    const MaskType& neighMask = DataAccess::get(permMask, neighPos);
                    // shifting signed integer is undefined behavior
                    const unsigned int nPermNeighbors = (unsigned int)neighMask >> 8;

                    // exclude non-periodic envelop
                    if((neighMask & PF_isPermeable) == false && nPermNeighbors > 0)
                    {
                        const T* sourceTensor = DataAccess::get(sourceTensorField, neighPos + offsetSource);

                        for(size_t iComp = 0; iComp < nComps; ++iComp)
                        {
                            const T& molesToDistribute = sourceTensor[iComp];
                            const T sourceVal = molesToDistribute / nPermNeighbors;

                            const Dot& offset = offsetVec[iComp];
                            Lattice& lattice = *static_cast<Lattice*>(fields[iComp+startLattices]);

                            plb::Cell<T,Descriptor>& cell = DataAccess::get(lattice, pos + offset);
                            T& source = *cell.getExternal(Descriptor<T>::ExternalField::scalarBeginsAt);

                            source += sourceVal;


                            const T& porosity = *cell.getExternal(Descriptor<T>::ExternalField::porosityBeginsAt);
                            auto concSlice = m_postReacConc[iComp];
                            U& conc = DataAccess::get(concSlice, absPos);

                            // std::cout << iComp << " conc: " << conc << " -> ";

                            const U term = sourceVal / porosity;
                            conc += term;

                            // std::cout << conc << std::endl;
                        }
                    }
                }
            }
        }
    }

private:
    const Dot& m_location;
    av::array_view<U, dim+1>& m_postReacConc;
};




}

#endif // SETNTENSORTOSOURCE_H
