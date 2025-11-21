#include "localaccessphreeqcrm.h"

#include "fielddecomposition.h"
#include <IPhreeqcPhast.h>
#include <CSelectedOutput.hxx>
#include <Solution.h>

#include "geometrytools.h"
#include "ippexception.h"
#include "scopedfloatingpointexception.h"

// debugging
#include "mpitools.h"
#include <chrono>
#include "ippstream.h"
#include "bench_tools.h"
//////

namespace IPP
{

LocalAccessPhreeqcRM::LocalAccessPhreeqcRM(int nxyz_arg,
                                           MPI_Comm parallelProcessData,
                                           PHRQ_io *io)
    : PhreeqcRM(nxyz_arg, parallelProcessData, io)
    , m_enabled(true)
    , m_decomp(nullptr)
    , m_isAllEnabled(false)
{

}

LocalAccessPhreeqcRM::~LocalAccessPhreeqcRM()
{
    m_decomp = nullptr;
}

void LocalAccessPhreeqcRM::init(const FieldDecomposition *decomp)
{
    MPI_CHECK_SYNC;

    assert(m_decomp == nullptr);
    m_decomp = decomp;
    applyDecomp(*decomp);

    m_saturationEnabled = saturation_worker;

    m_saturationCache.resize(m_saturationEnabled.size(), -1.0);
}

void LocalAccessPhreeqcRM::setGlobalEnabled(const bool state)
{
    m_enabled = state;
}

void LocalAccessPhreeqcRM::getEnabledCells(std::vector<char>& enabledCells) const
{
    enabledCells.resize(saturation_worker.size());
    for(size_t iCell = 0; iCell < m_saturationCache.size(); ++iCell)
    {
        if(saturation_worker[iCell] <= 0.0)
        {
            enabledCells[iCell] = false;
        }
        else
        {
            enabledCells[iCell] = true;
        }
    }
}

void LocalAccessPhreeqcRM::enableCell(const size_t iCellLocal)
{
    assert(m_decomp);
    IPPCheck::assertCheck(m_isAllEnabled == false);

    assert(saturation_worker.size() > iCellLocal);
    assert(m_saturationCache.size() == saturation_worker.size());

    assert(m_saturationCache[iCellLocal] > 0.0);
    std::swap(saturation_worker[iCellLocal], m_saturationCache[iCellLocal]);
    assert(m_saturationCache[iCellLocal] == -1.0);
    assert(saturation_worker[iCellLocal] > 0.0);
}

void LocalAccessPhreeqcRM::disableCell(const size_t iCellLocal)
{
    //    std::stringstream ss;
    //    ss << mpi_myself << " phreeqc disable: " << iCellLocal << std::endl;
    //    std::cout << ss.str();

    assert(m_decomp);
    IPPCheck::assertCheck(m_isAllEnabled == false);

    assert(saturation_worker.size() > iCellLocal);
    assert(m_saturationCache.size() == saturation_worker.size());

    assert(saturation_worker[iCellLocal] > 0.0);
    std::swap(m_saturationCache[iCellLocal], saturation_worker[iCellLocal]);
    assert(m_saturationCache[iCellLocal] > 0.0);
    assert(saturation_worker[iCellLocal] == -1.0);
}

IRM_RESULT LocalAccessPhreeqcRM::runCells()
{
    std::vector<int> r_vector;
    this->prepareErrorStreams(r_vector);

    this->runCellsLocal(r_vector);

    const IRM_RESULT status = this->retrieveErrors(r_vector);

    return status;
}

IRM_RESULT LocalAccessPhreeqcRM::runCellsIfEnabled()
{
    std::vector<int> r_vector;

    this->prepareErrorStreams(r_vector);


    if(m_enabled)
    {
        //        std::cout << MPIManager::getInstance().getRank() << " really run cells" << std::endl;
        this->runCellsLocal(r_vector);
    }
    else
    {
        //        std::cout << MPIManager::getInstance().getRank() << " omitting run cells" << std::endl;
    }

    MPI_CHECK_SYNC;

    const IRM_RESULT status = this->retrieveErrorsLocal(r_vector);

    MPI_CHECK_SYNC;

    return status;
}

IRM_RESULT LocalAccessPhreeqcRM::runCellsLocal()
{
    std::vector<int> r_vector;
    this->prepareErrorStreams(r_vector);

    this->runCellsLocal(r_vector);

    const IRM_RESULT status = this->retrieveErrorsLocal(r_vector);

    return status;
}

// static bool isValidPoros(const std::vector<double>& porosities)
// {
//     for(const double& poros : porosities)
//     {
//         if(poros < 0.0 || poros > 1.0)
//         {
//             return false;
//         }
//     }

//     return true;
// }

void LocalAccessPhreeqcRM::setPorosityLocal(const std::vector<double>& porosities)
{
    assert(porosity_worker.size() == porosities.size());

    //    std::stringstream ss;
    //    ss << mpi_myself << " poros: " << std::endl;
    //    for(size_t i = 0; i < porosities.size(); ++i)
    //    {
    //        ss << "\t" << i << "\t" << porosities[i] << std::endl;
    //    }
    //    std::cout << ss.str();


    porosity_worker = porosities;
}

const std::vector<double>& LocalAccessPhreeqcRM::getPorosityLocal() const
{
    return porosity_worker;
}

void LocalAccessPhreeqcRM::getSolutionVolumeLocal(std::vector<double>& volume) const
{
    const size_t nCellsLocal = m_decomp->getLocalNumberOfCells();

    volume.resize(nCellsLocal, INACTIVE_CELL_VALUE);


    int n = this->mpi_myself;
    const size_t startCell = this->start_cell[n];
    assert(this->end_cell[n] >= (int) startCell);
    assert((this->end_cell[n] - startCell + 1) == nCellsLocal);

    IPhreeqcPhast &internalPhreeqc = *workers[0];

    // fill solution_volume
    for (size_t iCellLocal = 0; iCellLocal < nCellsLocal; ++iCellLocal)
    {
        const size_t iCellGlobal = startCell + iCellLocal;

        const cxxSolution* solution = internalPhreeqc.Get_solution(iCellGlobal);
        volume[iCellLocal] = solution->Get_soln_vol();
    }
}

void LocalAccessPhreeqcRM::enableAllCellsTemp()
{
    assert(m_decomp);
    IPPCheck::assertCheck(m_isAllEnabled == false);

    assert(m_saturationEnabled.size() == saturation_worker.size());
    saturation_worker.swap(m_saturationEnabled);
    m_isAllEnabled = true;
}

void LocalAccessPhreeqcRM::resetToEnabledStatus()
{
    assert(m_decomp);
    IPPCheck::assertCheck(m_isAllEnabled);

    assert(m_saturationEnabled.size() == saturation_worker.size());
    saturation_worker.swap(m_saturationEnabled);
    m_isAllEnabled = false;
}

void LocalAccessPhreeqcRM::prepareErrorStreams(std::vector<int>& r_vector)
{
    IPhreeqcPhast * phast_iphreeqc_worker = this->workers[0];
    // phast_iphreeqc_worker->Get_PhreeqcPtr()->Set_run_cells_one_step(true);


    delete phast_iphreeqc_worker->Get_out_stream();
    phast_iphreeqc_worker->Set_out_stream(new std::ostringstream);

    delete phast_iphreeqc_worker->Get_punch_stream();
    phast_iphreeqc_worker->Set_punch_stream(new std::ostringstream);


    r_vector.resize(1);
}

void LocalAccessPhreeqcRM::runCellsLocal(std::vector<int>& r_vector)
{
    this->phreeqcrm_error_string.clear();



    // check that all solutions are defined
    //    if (this->need_error_check)
    //    {
    //        this->need_error_check = false;
    //        try
    //        {
    //            CheckCells();
    //        }
    //        catch(...)
    //        {
    //            MPI_Barrier(phreeqcrm_comm);
    //            return this->ReturnHandler(IRM_FAIL, "PhreeqcRM::RunCells");
    //        }
    //    }



    {
        // for some reason phreeqc sometimes produces FPE INV
        // additioanlly there is a div zero quite frequently in model.cpp:3400-ish
        // latter one was tried to fix with tolerance check within phreeqc itself
        // but was disabled again, due to possible unexpected side effects

        ScopedDisableFloatingPointException fpe;
        r_vector[0] = RunCellsThread(0);
    }

    if (this->partition_uz_solids)
    {
        old_saturation_worker = saturation_worker;
        if (mpi_myself == 0)
        {
            old_saturation_root = saturation_root;
        }
    }

}


IRM_RESULT LocalAccessPhreeqcRM::retrieveErrorsLocal(std::vector<int>& r_vector)
{
    IRM_RESULT return_value = IRM_OK;


    // Count errors and write error messages
    try
    {
        // Check for errors
        this->error_count = 0;

        // Write error messages
        for (size_t n = 0; n < r_vector.size(); n++)
        {
            if (r_vector[n] != 0)
            {
                // print error
                std::ostringstream e_stream;
                e_stream << "Process " << mpi_myself << std::endl;
                this->ErrorMessage(e_stream.str());
                this->ErrorMessage(this->workers[n]->GetErrorString(), false);
                this->error_count++;
            }
        }
        if (error_count > 0)
            throw PhreeqcRMStop();

        return IRM_OK;
    }
    catch (...)
    {
        return_value = IRM_FAIL;
    }

    return this->ReturnHandler(return_value, "PhreeqcRM::RunCells");
}

IRM_RESULT LocalAccessPhreeqcRM::retrieveErrors(std::vector<int>& r_vector)
{
    IRM_RESULT return_value = IRM_OK;
    std::vector<char> char_buffer;

    // write output results
    if (this->print_chemistry_on[0])
    {
        for (int n = 0; n < this->mpi_tasks; n++)
        {
            // Need to transfer output stream to root and print
            if (this->mpi_myself == n)
            {
                if (n == 0)
                {
                    IPhreeqcPhast* worker = this->workers[0];
                    std::ostringstream* stream = worker->Get_out_stream();
                    this->OutputMessage(stream->str().c_str());
                    delete stream;
                    worker->Set_out_stream(NULL);
                }
                else
                {

                    int size = (int) this->workers[0]->Get_out_stream()->str().size();
                    MPI_Send(&size, 1, MPI_INT, 0, 0, phreeqcrm_comm);

                    MPI_Send((void *) this->workers[0]->Get_out_stream()->str().c_str(), size, MPI_CHAR, 0, 0, phreeqcrm_comm);
                    delete this->workers[0]->Get_out_stream();
                    this->workers[0]->Set_out_stream(NULL);
                }
            }
            else if (this->mpi_myself == 0)
            {
                MPI_Status mpi_status;
                int size;
                MPI_Recv(&size, 1, MPI_INT, n, 0, phreeqcrm_comm, &mpi_status);
                char_buffer.resize(size + 1);
                MPI_Recv((void *) &char_buffer.front(), size, MPI_CHAR, n, 0, phreeqcrm_comm, &mpi_status);
                char_buffer[size] = '\0';
                this->OutputMessage(&char_buffer.front());
            }
        }
    }

    // Count errors and write error messages
    try
    {
        HandleErrorsInternal(r_vector);

#ifdef IPP_DEBUG
        this->CheckSelectedOutput();
#endif
        // do not touch rebalancing

        // Rebalance load
        //        double t0 = (double) MPI_Wtime();
        //        this->RebalanceLoad();

    }
    catch (...)
    {
        return_value = IRM_FAIL;
    }

    return this->ReturnHandler(return_value, "PhreeqcRM::RunCells");
}

IRM_RESULT LocalAccessPhreeqcRM::getConcentrationsLocal(std::vector<double> &c)
{

    assert(m_decomp);

    this->phreeqcrm_error_string.clear();
    IRM_RESULT return_value = IRM_OK;

    // convert Reaction module solution data to concentrations for transport

    const int n = mpi_myself;
    const size_t nCells = (1+end_cell[n]) - start_cell[n];
    const size_t nComps = components.size();
    // std::cout << nCells*nComps << std::endl;
    m_solns.resize(nCells * nComps);

    try
    {
        size_t index = 0;

        // Put solutions into a vector

        for (int iCell = start_cell[n]; iCell <= end_cell[n]; iCell++)
        {
            // load fractions into d
            cxxSolution * cxxsoln_ptr = this->GetWorkers()[0]->Get_solution(iCell);
            assert (cxxsoln_ptr);

            double v, dens;

            if (this->use_solution_density_volume)
            {
                v = cxxsoln_ptr->Get_soln_vol();
                dens = cxxsoln_ptr->Get_density();
            }
            else
            {
                //int k = this->backward_mapping[j][0];
                int l = iCell - start_cell[n];
                v = saturation_worker[l] * porosity_worker[l] * rv_worker[l];

                if (v <= 0)
                {
                    v = cxxsoln_ptr->Get_soln_vol();
                }
                dens = density_worker[l];
            }

            std::vector<double> d;
            this->cxxSolution2concentration(cxxsoln_ptr, d, v, dens);

            for (size_t iComp = 0; iComp < components.size(); iComp++)
            {
                const double& conc = d[iComp];
                m_solns[index] = conc;
                ++index;
            }
        }



        GeometryTools::transpose(m_solns, nComps, nCells);

        c.swap(m_solns);

    }
    catch (...)
    {
        return_value = IRM_FAIL;
    }
    return this->ReturnHandler(return_value, "PhreeqcRM::GetConcentrations");

}

IRM_RESULT LocalAccessPhreeqcRM::setConcentrationsLocal(const std::vector<double> &c)
{
    assert(m_decomp);

    // stupid PhreeqcRM interface needs full set / global of concentrations
    // but Concentrations2Solutions function uses only the relevant part
    // OLD: initializing the not used cells with INACTIVE_CELL_VALUE (hopefully this crashes upon error...)
    // REMARK: values are now uninitialized due to performance benefits. scratch space is allocated only once
    // WARNING: end is NOT end in c++ terminology. in fact end is last in phreeqc

    const size_t beginCell = start_cell[mpi_myself];
    const size_t nCellsLocal = (1+end_cell[mpi_myself]) - beginCell;
    const size_t nComps = components.size();
    const size_t dstIndexOff = beginCell * nComps;

    c_chem.resize(nxyz * nComps, INACTIVE_CELL_VALUE);

    for (size_t iCellLocal = 0; iCellLocal < nCellsLocal; ++iCellLocal)
    {
        // TODO: disable backward mapping?
        //int j = this->backward_mapping[iCellLocal][0];

        const size_t compOff = dstIndexOff + nComps * iCellLocal;

        for (size_t iComp = 0; iComp < nComps; iComp++)
        {
            const size_t srcIndex = iComp * nCellsLocal + iCellLocal;
            const size_t dstIndex = compOff + iComp;
            const double& srcVal = c[srcIndex];

            // slightly negative concentrations are allowed
            // limitation in phreeqcrm was lifted
            // assert(iComp == 3 || srcVal >= 0.0);

            c_chem[dstIndex] = srcVal;
        }
    }

    // saturation must be enabled temporarily otherwise the inactive cells will not be updated
    this->enableAllCellsTemp();

    {
        // stupid Concentrations2Solutions function does trigger fpe....
        ScopedDisableFloatingPointException fpe;
        const int worker = 0;
        this->Concentrations2Solutions(worker, c_chem);
    }

    this->resetToEnabledStatus();

    IRM_RESULT return_value = IRM_OK;
    return this->ReturnHandler(return_value, "PhreeqcRM::SetConcentrations");
}

IRM_RESULT LocalAccessPhreeqcRM::getSelectedOutputLocal(std::vector<double> &so)
{

    assert(m_decomp);

    this->phreeqcrm_error_string.clear();
    IRM_RESULT return_value = IRM_OK;

    try
    {
        const IPhreeqcPhast* worker = workers[0];
        int n_user = worker->GetCurrentSelectedOutputUserNumber();

        // no bcast for local operations!
        // MPI_Bcast(&n_user,  1, MPI_INT, 0, phreeqcrm_comm);
        // if (n_user < 0)
        // {
        //     this->ErrorHandler(IRM_INVALIDARG, "No selected output defined");
        // }

        std::map<int, CSelectedOutput>::const_iterator it = worker->CSelectedOutputMap.find(n_user);
        if (it == worker->CSelectedOutputMap.end() || this->SetCurrentSelectedOutputUserNumber(n_user) < 0)
        {
            this->ErrorHandler(IRM_INVALIDARG, "Selected output not found");
        }

        const CSelectedOutput& cso = it->second;

        int ncol = this->GetSelectedOutputColumnCount();


        // fill with INACTIVE_CELL_VALUE

        int nrowLocal;
        std::vector<double> dbuffer;
        cso.Doublize(nrowLocal, ncol, dbuffer);

        assert(dbuffer.size() == m_decomp->getLocalNumberOfCells() * ncol);

        // TODO: check for backward mapping
        so.swap(dbuffer);

    }
    catch (...)
    {
        return_value = IRM_FAIL;
    }
    return this->ReturnHandler(return_value, "PhreeqcRM::GetSelectedOutput");
}

void LocalAccessPhreeqcRM::setTimeLocal(const double &t)
{
    time = t;
}

void LocalAccessPhreeqcRM::setFilePrefixLocal(const std::string &prefix)
{
    assert(prefix.empty() == false);
    this->file_prefix = prefix;
}

void LocalAccessPhreeqcRM::setPrintChemistryMaskLocal(std::vector<int> &m)
{
    this->print_chem_mask_root = m;
}

void LocalAccessPhreeqcRM::setPrintChemistryOnLocal(bool worker, bool ip, bool utility)
{
    this->print_chemistry_on[0] = worker;
    this->print_chemistry_on[1] = ip;
    this->print_chemistry_on[2] = utility;
}



static void collectToRun(const int nthreads,
                         bool runInitialPhreeqc,
                         bool runWorkers,
                         bool runUtility,
                         std::vector<bool>& run)
{
    run.resize(nthreads + 2, false);

    if(runWorkers)
    {
        for (int i = 0; i < nthreads; i++)
        {
            run[i] = true;
        }
    }
    if(runInitialPhreeqc)
    {
        run[nthreads] = true;
    }
    if (runUtility)
    {
        run[nthreads + 1] = true;
    }
}

IRM_RESULT LocalAccessPhreeqcRM::collectErrors(const std::vector<IRM_RESULT>& errors)
{
    IRM_RESULT finalResult = IRM_OK;
    for (int n = 0; n < nthreads + 2; n++)
    {
        if(errors[n] != IRM_OK)
        {
            finalResult = errors[n];
            IPhreeqcPhast* worker = workers[n];
            const char* err = worker->GetErrorString();
            std::cout << err << std::endl;
        }
    }

    return finalResult;
}

IRM_RESULT LocalAccessPhreeqcRM::runStringLocal(bool runWorkers,
                                                bool runInitialPhreeqc,
                                                bool runUtility,
                                                const std::string &input)
{
    std::vector<bool> run;
    collectToRun(nthreads, runInitialPhreeqc, runWorkers, runUtility, run);

    std::vector<IRM_RESULT> errors(nthreads + 2, IRM_OK);
    for (int n = 0; n < nthreads + 2; n++)
    {
        if (run[n])
        {
            errors[n] = this->RunStringThread(n, input);
        }
    }

    const IRM_RESULT finalResult = collectErrors(errors);

    return finalResult;
}

IRM_RESULT LocalAccessPhreeqcRM::runStringLocal(bool runWorkers,
                                                bool runInitialPhreeqc,
                                                bool runUtility,
                                                std::istream& input)
{
    std::vector<bool> run;
    collectToRun(nthreads, runInitialPhreeqc, runWorkers, runUtility, run);

    std::vector<IRM_RESULT> errors(nthreads + 2, IRM_OK);
    for (int n = 0; n < nthreads + 2; n++)
    {
        if (run[n])
        {
            errors[n] = this->runStringStreamThread(n, input);
        }
    }

    const IRM_RESULT finalResult = collectErrors(errors);

    return finalResult;
}

IRM_RESULT LocalAccessPhreeqcRM::runStringStreamThread(int n, std::istream& input)

{
    try
    {
        IPhreeqcPhast * iphreeqc_phast_worker = this->GetWorkers()[n];

        iphreeqc_phast_worker->SetOutputFileOn(false);
        iphreeqc_phast_worker->SetErrorFileOn(false);
        iphreeqc_phast_worker->SetLogFileOn(false);
        iphreeqc_phast_worker->SetSelectedOutputStringOn(false);
        iphreeqc_phast_worker->SetSelectedOutputFileOn(false);

        // Set output string on
        if (n < this->nthreads)
        {
            iphreeqc_phast_worker->SetOutputStringOn(this->print_chemistry_on[0]);
        }
        else if (n == this->nthreads)
        {
            iphreeqc_phast_worker->SetOutputStringOn(this->print_chemistry_on[1]);
        }
        else
        {
            iphreeqc_phast_worker->SetOutputStringOn(this->print_chemistry_on[2]);
        }
        // Run chemistry file
        if (iphreeqc_phast_worker->RunStringStream(input) > 0)
        {
            this->ErrorMessage(iphreeqc_phast_worker->GetErrorString());
            throw PhreeqcRMStop();
        }

        if (iphreeqc_phast_worker->GetOutputStringOn())
        {
            this->OutputMessage(iphreeqc_phast_worker->GetOutputString());
        }
    }
    catch (PhreeqcRMStop& e)
    {
        return IRM_FAIL;
    }
    catch (...)
    {
        std::ostringstream e_stream;
        e_stream << "RunString failed in worker " << n << "from an unhandled exception.\n";
        this->ErrorMessage(e_stream.str());
        return IRM_FAIL;
    }
    return IRM_OK;
}

size_t LocalAccessPhreeqcRM::getNumberCells() const
{
    const size_t result = (1 + end_cell.at(mpi_myself)) - start_cell.at(mpi_myself);
    return result;
}

void LocalAccessPhreeqcRM::applyDecomp(const FieldDecomposition &decomp)
{
    const std::vector<IPPBox3DLong>& domains = decomp.getDomains();
    assert((int)domains.size() == mpi_tasks);

    std::vector<int> nCellsVec;
    for(size_t iProc = 0; iProc < domains.size(); ++iProc)
    {
        const IPPBox3DLong& domain = domains[iProc];
        const IPPVector3DLong size = domain.getSize();
        const int nCells = size[0] * size[1] * size[2];
        nCellsVec.push_back(nCells);
    }

    std::vector<int>& new_start_cell = start_cell;
    std::vector<int>& new_end_cell = end_cell;
    new_start_cell.clear();
    new_end_cell.clear();


    int cellBegin = 0;

    for (int i = 0; i < mpi_tasks; i++)
    {
        const int nCells = nCellsVec[i];

        const int cellEnd = cellBegin + nCells;
        const int cellLast = cellEnd - 1;

        new_start_cell.push_back(cellBegin);
        new_end_cell.push_back(cellLast);

        cellBegin += nCells;
    }


    ScatterNchem(print_chem_mask_root, print_chem_mask_worker);
    ScatterNchem(density_root, density_worker);
    ScatterNchem(porosity_root, porosity_worker);
    ScatterNchem(rv_root, rv_worker);
    ScatterNchem(saturation_root, saturation_worker);

}

}
