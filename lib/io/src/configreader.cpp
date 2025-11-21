#include "configreader.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>

#include <json/json.h>

#include <arraydimensionconvert.h>
#include <bcgenerator.h>
#include <configdatawrapper.h>
#include <ippexception.h>
#include <phasenametodissprecipbehaviour.h>
#include <scopedfloatingpointexception.h>

#include "cemhyd3dimagereader.h"
#include "geometrytools.h"
#include "initfunctors.h"
#include "inputgeometryscale.h"
#include "saturationindexcalcfactory.h"
#include "vtkimagereader.h"

#include "abstractmultiscalediffusioncalc.h"
#include "abstractporositycalc.h"
#include "abstractsicalc.h"

#include "dissolveonlyonbelowthreshold.h"
#include "nucleationdissolveonly.h"

// dummys
#include "nonreturningbcgen.h"
#include "nonreturningmultiscalediffusionfactory.h"
#include "nonreturningporositymodelfactory.h"
#include "nonreturningsicalcfactory.h"

namespace IPP {

static ResultType getResultType(const std::string &name)
{
    if (boost::iequals(name, "porosity")) {
        return RT_Porosity;
    } else if (boost::iequals(name, "capillary_porosity")) {
        return RT_CapillaryPorosity;
    } else if (boost::iequals(name, "inert_frac")) {
        return RT_InertFrac;
    } else if (boost::iequals(name, "diff_coef")) {
        return RT_DiffusionCoef;
    } else if (boost::iequals(name, "diff_transport_scalar")) {
        return RT_DiffusiveTransportScalar;
    } else if (boost::iequals(name, "dist_field")) {
        return RT_DistanceField;
    } else if (boost::iequals(name, "enabled_cells")) {
        return RT_EnabledCells;
    } else if (boost::iequals(name, "phases")) {
        return RT_Phase;
    } else if (boost::iequals(name, "LBM_velocity")) {
        return RT_LBM_Velocity;
    } else if (boost::iequals(name, "LBM_gradient")) {
        return RT_LBM_Gradient;
    } else if (boost::iequals(name, "LBM_post_transport_conc")) {
        return RT_LBM_PostTransportConc;
    } else if (boost::iequals(name, "LBM_source")) {
        return RT_LBM_Source;
    } else if (boost::iequals(name, "LBM_trans_poros")) {
        return RT_LBM_TransPoros;
    } else if (boost::iequals(name, "Budget")) {
        return RT_Budget;
    } else if (boost::iequals(name, "perm_info")) {
        return RT_PermeabilityInfo;
    } else {
        throw std::runtime_error("unknown result type: " + name);
    }
}

void parseResults(const Json::Value &root, IPPConfig &config)
{
    IPPConfig::Results &results = config.results;
    ResultsToWrite &resultsToWrite = results.resultsToWrite;

    // adding default results
    resultsToWrite.add(RT_ComponentConcentration);
    // resultsToWrite.add(RT_ComponentFlux);

    if (config.isChemistryEnabled) {
        resultsToWrite.add(RT_Phase);
    }

    if (root.isMember("results")) {
        const Json::Value resultsNode = root["results"];

        const Json::Value enabledArrNode = resultsNode["results_arrays_enabled"];
        for (size_t i = 0; i < enabledArrNode.size(); ++i) {
            const Json::Value enabledNode = enabledArrNode[(int) i];
            const std::string name = enabledNode.asString();

            const ResultType result = getResultType(name);
            resultsToWrite.add(result);
        }

        const Json::Value disabledArrNode = resultsNode["results_arrays_disabled"];
        for (size_t i = 0; i < disabledArrNode.size(); ++i) {
            const Json::Value disabledNode = disabledArrNode[(int) i];
            const std::string name = disabledNode.asString();

            const ResultType result = getResultType(name);
            resultsToWrite.remove(result);
        }

        // initialize defaults from flux_convergence_tolerance
        const double uniqueConvTol = resultsNode.get("flux_convergence_tolerance", 0.0).asDouble();
        config.fluxConvergenceTolLow = uniqueConvTol;
        config.fluxConvergenceTolHigh = uniqueConvTol;

        if (resultsNode.isMember("flux_convergence_tolerance_low")) {
            config.fluxConvergenceTolLow = resultsNode["flux_convergence_tolerance_low"].asDouble();
        }

        if (resultsNode.isMember("flux_convergence_tolerance_high")) {
            config.fluxConvergenceTolHigh = resultsNode["flux_convergence_tolerance_high"].asDouble();
        }

        config.fluxConvergenceTolFactor = resultsNode.get("flux_convergence_tolerance_factor", 2.0)
                                              .asDouble();

        const Json::Value &diffEffNode = resultsNode["diff_eff"];
        for (int iDim = 0; iDim < (int) diffEffNode.size(); ++iDim) {
            const Json::Value &node = diffEffNode[iDim];

            const size_t dim = node["dim"].asUInt();
            const std::string tracerName = node["tracer_name"].asString();

            if (std::find_if(results.diffEffInfos.begin(),
                             results.diffEffInfos.end(),
                             [&](const IPPConfig::Results::DiffEffInfo &info) {
                                 return info.dim == dim;
                             })
                    == results.diffEffInfos.end()

                &&

                std::find_if(results.diffEffInfos.begin(),
                             results.diffEffInfos.end(),
                             [&](const IPPConfig::Results::DiffEffInfo &info) {
                                 return info.tracerName == tracerName;
                             })
                    == results.diffEffInfos.end()) {
                results.diffEffInfos.push_back(IPPConfig::Results::DiffEffInfo(dim, tracerName));
            } else {
                throw std::runtime_error(
                    "inert tracer for effective diffusion defined more than once: "
                    + std::to_string(dim) + " " + tracerName);
            }
        }

        const Json::Value &auxNode = resultsNode["aux_data"];
        for (size_t i = 0; i < auxNode.size(); ++i) {
            const Json::Value &entry = auxNode[(int) i];
            const std::vector<std::string> &names = entry.getMemberNames();
            const std::string name = names.front();
            const std::string command = entry[name].asString();
            results.auxData.push_back(IPPConfig::Results::NameCommand(name, command));
        }
    }
}

// input geom is invalidated
void applyScaling(const size_t nxSrc,
                  const size_t nySrc,
                  const size_t nzSrc,
                  const size_t scalePasses,
                  std::vector<size_t> &inputGeom,
                  std::vector<size_t> &output)
{
    InputGeometryScale::ColumnMajorArrayDimensionConvert srcIndexConv(nxSrc, nySrc, nzSrc);

    for (size_t iPass = 0; iPass < scalePasses; ++iPass) {
        const size_t nxDst = InputGeometryScale::calcUpscaledSize(srcIndexConv.getNx(), 1);
        const size_t nyDst = InputGeometryScale::calcUpscaledSize(srcIndexConv.getNy(), 1);
        const size_t nzDst = InputGeometryScale::calcUpscaledSize(srcIndexConv.getNz(), 1);
        InputGeometryScale::ColumnMajorArrayDimensionConvert dstIndexConv(nxDst, nyDst, nzDst);

        std::vector<size_t> scaledResult;
        InputGeometryScale::scale(srcIndexConv, dstIndexConv, inputGeom, scaledResult);

        srcIndexConv = dstIndexConv;
        inputGeom.swap(scaledResult);
    }

    output.swap(inputGeom);
}

typedef std::unordered_map<size_t, Composition *> IDtoComposition;
void setupGeometry(const size_t nxSrc,
                   const size_t nySrc,
                   const size_t nzSrc,
                   const IDtoComposition &idToComposition,
                   const IPPVector3DInt bounds,
                   const IPPVector3DInt origin,
                   const size_t scalePasses,
                   const ArrayDimensionConvert inputIndexConv,
                   const ArrayDimensionConvert outputindexConv,
                   std::vector<size_t> &inputGeom,
                   IPPConfig &config)
{
    // apply scaling
    std::vector<size_t> domainsVec;
    applyScaling(nxSrc, nySrc, nzSrc, scalePasses, inputGeom, domainsVec);

    std::vector<size_t> uniqueDomains = domainsVec;
    {
        std::sort(uniqueDomains.begin(), uniqueDomains.end());
        std::vector<size_t>::const_iterator last = std::unique(uniqueDomains.begin(),
                                                               uniqueDomains.end());
        uniqueDomains.erase(last, uniqueDomains.end());
    }

    std::vector<Domain> &domains = config.domains;

    assert(domains.empty());
    std::unordered_map<size_t, size_t> domainIDtoIndex;
    for (size_t i = 0; i < uniqueDomains.size(); ++i) {
        const int domainID = uniqueDomains[i];

        if (idToComposition.find(domainID) == idToComposition.end()) {
            throw std::runtime_error("domain with id not defined in ID_Mapping: "
                                     + std::to_string(domainID));
        }

        const Composition *compPtr = idToComposition.at(domainID);
        domains.push_back(Domain(compPtr));
        domainIDtoIndex[domainID] = i;
    }

    for (size_t iCell = 0; iCell < domainsVec.size(); ++iCell) {
        const size_t domainID = domainsVec[iCell];
        const size_t index = domainIDtoIndex[domainID];

        const IPPVector3DInt pos = inputIndexConv.colMajorToRowMajorCoord(iCell);
        const IPPVector3DInt absPos = pos + origin;

        // just crop if its outside of domain
        if (GeometryTools::isWithin(absPos, bounds)) {
            const size_t absIndex = outputindexConv.calcIndex(absPos);
            Domain &domainConf = domains[index];
            domainConf.cells.push_back(absIndex);
        }
    }

    for (size_t iDomain = 0; iDomain < domains.size(); ++iDomain) {
        Domain &domain = domains[iDomain];
        std::vector<int> &cells = domain.cells;
        const size_t oldSize = cells.size();

        std::sort(cells.begin(), cells.end());
        auto last = std::unique(cells.begin(), cells.end());
        cells.erase(last, cells.end());

        IPPCheck::assertCheck(cells.size() == oldSize, "cells do not fit");
    }
}

void readGeometry(const Json::Value geometry,
                  const IDtoComposition &idToComposition,
                  IPPConfig &config)
{
    const IPPVector3DInt bounds = {{(int) config.nx, (int) config.ny, (int) config.nz}};

    const Json::Value &dims = geometry["dims"];
    const size_t arrNx = dims.get(0u, config.nx).asUInt();
    const size_t arrNy = dims.get(1u, config.ny).asUInt();
    const size_t arrNz = dims.get(2u, config.nz).asUInt();

    const Json::Value &originNode = geometry["origin"];
    const IPPVector3DInt origin = {{originNode.get(0u, 0).asInt(),
                                    originNode.get(1u, 0).asInt(),
                                    originNode.get(2u, 0).asInt()}};

    const ArrayDimensionConvert inputIndexConv(arrNx, arrNy, arrNz);
    const ArrayDimensionConvert outputIndexConv(bounds[0], bounds[1], bounds[2]);

    // reading geometric informations

    const size_t scalePasses = geometry.get("scale_passes", 0).asUInt64();

    const size_t nxSrc = InputGeometryScale::calcDownscaledSize(arrNx, scalePasses);
    const size_t nxUpscaled = InputGeometryScale::calcUpscaledSize(nxSrc, scalePasses);
    IPPCheck::assertCheck(arrNx == nxUpscaled,
                          "nx is not compatible to scale passes: " + std::to_string(arrNx) + " vs. "
                              + std::to_string(nxUpscaled));

    const size_t nySrc = InputGeometryScale::calcDownscaledSize(arrNy, scalePasses);
    const size_t nyUpscaled = InputGeometryScale::calcUpscaledSize(nySrc, scalePasses);
    IPPCheck::assertCheck(arrNy == nyUpscaled,
                          "ny is not compatible to scale passes: " + std::to_string(arrNy) + " vs. "
                              + std::to_string(nyUpscaled));

    const size_t nzSrc = InputGeometryScale::calcDownscaledSize(arrNz, scalePasses);
    const size_t nzUpscaled = InputGeometryScale::calcUpscaledSize(nzSrc, scalePasses);
    IPPCheck::assertCheck(arrNz == nzUpscaled,
                          "nz is not compatible to scale passes: " + std::to_string(arrNz) + " vs. "
                              + std::to_string(nzUpscaled));

    const Json::Value arrayNode = geometry["array"];

    if (arrayNode.isArray()) {
        std::vector<size_t> inputGeom(arrayNode.size());
        for (size_t iCell = 0; iCell < arrayNode.size(); ++iCell) {
            const size_t domainID = arrayNode[(int) iCell].asUInt64();
            inputGeom[iCell] = domainID;
        }

        setupGeometry(nxSrc,
                      nySrc,
                      nzSrc,
                      idToComposition,
                      bounds,
                      origin,
                      scalePasses,
                      inputIndexConv,
                      outputIndexConv,
                      inputGeom,
                      config);

    } else if (geometry.isMember("cemhyd_file")) {
        const std::string cemhydFile = geometry["cemhyd_file"].asString();

        std::vector<size_t> input;
        Cemhyd3DImageReader::read(cemhydFile, input);

        std::vector<size_t> scaled;

        if (scalePasses > 0) {
            // fill to nXYZ;
            const size_t nxyz = nxSrc * nySrc * nzSrc;
            for (size_t i = input.size(); i < nxyz; ++i) {
                input.push_back(0);
            }

            applyScaling(nxSrc, nySrc, nzSrc, scalePasses, input, scaled);
        } else {
            scaled.swap(input);
        }

        Cemhyd3DImageReader::convert(inputIndexConv,
                                     idToComposition,
                                     scaled,
                                     origin,
                                     bounds,
                                     config.domains);

    } else if (geometry.isMember("vtk_file")) {
        const std::string vtkFile = geometry["vtk_file"].asString();

        std::vector<size_t> inputGeom;
        VTKImageReader::read(vtkFile, inputGeom);

        setupGeometry(nxSrc,
                      nySrc,
                      nzSrc,
                      idToComposition,
                      bounds,
                      origin,
                      scalePasses,
                      inputIndexConv,
                      outputIndexConv,
                      inputGeom,
                      config);
    } else {
        throw std::runtime_error("no geometry array defined");
    }

    // fix undefined cells
    if (geometry.isMember("undefined_cells_domain")) {
        const size_t nxyz = config.nx * config.ny * config.nz;
        std::vector<char> isDefined(nxyz, false);

        for (const Domain &dom : config.domains) {
            for (size_t iCell : dom.cells) {
                assert(isDefined.size() > iCell);
                isDefined[iCell] = true;
            }
        }

        const std::vector<char>::const_iterator it = std::find(isDefined.begin(),
                                                               isDefined.end(),
                                                               false);
        if (it != isDefined.end()) {
            // at least one cell is not defined, adding new domain

            const size_t compId = geometry["undefined_cells_domain"].asUInt();
            const Composition *comp = idToComposition.at(compId);
            config.domains.push_back(Domain(comp));
            Domain &dom = config.domains.back();

            for (size_t iCell = 0; iCell < isDefined.size(); ++iCell) {
                if (isDefined[iCell] == false) {
                    dom.cells.push_back(iCell);
                }
            }
        }
    }
}

static LatticeBoltzmannCollisionType getCollisionType(const Json::Value &root)
{
    const std::string collTypeStr = root.get("LBCollType", "PTRT").asString();
    if (collTypeStr == "DVSRT") {
        return LBCT_DVSRT;
    } else if (collTypeStr == "PTRT") {
        return LBCT_PTRT;
    } else {
        throw std::runtime_error("unknown collision type: " + collTypeStr);
    }
}

static PorosityThresholds parsePorosityThresholds(const bool is3DSimulation,
                                                  const double &porosMinThresh,
                                                  const Json::Value &doNode)
{
    const double porosUpperThresh = doNode.get("poros_upper_thresh", 1.0).asDouble();

    const double edgeThresh = 0.7;
    const double cornerThresh = is3DSimulation ? 1.0 - 0.5 / 6.0
                                               : -1.0; // there is no "corner" in 2D

    const PorosityThresholds thresh = {porosMinThresh, porosUpperThresh, edgeThresh, cornerThresh};
    return thresh;
}

static void parseNucleationData(const Json::Value &doNode, NucleationData &nameToData)
{
    IPPCheck::assertCheck(doNode.isMember("nucleation_infos"));
    const Json::Value &nucleation_infos = doNode["nucleation_infos"];

    for (size_t iNode = 0; iNode < nucleation_infos.size(); ++iNode) {
        const Json::Value &phaseNode = nucleation_infos[(int) iNode];

        IPPCheck::assertCheck(phaseNode.isMember("phase"));
        const std::string phaseName = phaseNode.get("phase", "").asString();

        if (nameToData.find(phaseName) != nameToData.end()) {
            throw std::runtime_error("nucleation data for phase defined multiple times: "
                                     + phaseName);
        }
        NucleationDissolveOnlyData &data = nameToData[phaseName];

        IPPCheck::assertCheck(phaseNode.isMember("monomer"));
        data.monomer = phaseNode["monomer"].asString();

        IPPCheck::assertCheck(phaseNode.isMember("molecule_vol"));
        data.V_molecule = phaseNode["molecule_vol"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("D"));
        data.D = phaseNode["D"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("T"));
        data.T = phaseNode["T"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("HON_IT_water"));
        data.HON_IT_water = phaseNode["HON_IT_water"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("HEN_IT_substrate"));
        data.HEN_IT_substrate = phaseNode["HEN_IT_substrate"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("HEN_IT_self"));
        data.HEN_IT_self = phaseNode["HEN_IT_self"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("HON_N0_water"));
        data.HON_N0_water = phaseNode["HON_N0_water"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("HEN_N0_substrate"));
        data.HEN_N0_substrate = phaseNode["HEN_N0_substrate"].asDouble();

        IPPCheck::assertCheck(phaseNode.isMember("HEN_N0_self"));
        data.HEN_N0_self = phaseNode["HEN_N0_self"].asDouble();
    }
}

static AbstractDissolveOnlyCalc *parseDissolveOnlyFunc(const bool is3DSimulation,
                                                       const Json::Value &doNode,
                                                       const double &porosThreshDefault)
{
    const std::string technique = doNode["technique"].asString();
    const double porosMinThresh = doNode.get("poros_lower_thresh", porosThreshDefault).asDouble();

    if (technique == "min_poros") {
        return new BelowThreshDissolveOnly(porosMinThresh);
    } else if (technique == "neigh_aware") {
        const PorosityThresholds thresh = parsePorosityThresholds(is3DSimulation,
                                                                  porosMinThresh,
                                                                  doNode);
        return new BelowThreshNeighborAwareDissolveOnly(thresh);
    } else if (technique == "nucleation") {
        const PorosityThresholds thresh = parsePorosityThresholds(is3DSimulation,
                                                                  porosMinThresh,
                                                                  doNode);

        NucleationDissolveOnly *nuclCalc = new NucleationDissolveOnly(thresh);
        NucleationData &nameToData = nuclCalc->getNucleationData();

        parseNucleationData(doNode, nameToData);

        return nuclCalc;
    } else {
        throw std::runtime_error("unknown dissolve only technique: " + technique);
    }
}

static IPPBox3DInt readBox(const Json::Value &node, const size_t dims)
{
    IPPCheck::assertCheck(node.size() == dims * 2, "define box as: [ x0, x1, y0, y1 (,z0, z2) ]");

    if (dims == 2) {
        const int x0 = node[0].asInt();
        const int x1 = node[1].asInt();
        const int y0 = node[2].asInt();
        const int y1 = node[3].asInt();
        const IPPVector3DInt lower = {{x0, y0, 0}};
        const IPPVector3DInt upper = {{x1, y1, 0}};
        return IPPBox3DInt(lower, upper);
    } else {
        const int x0 = node[0].asInt();
        const int x1 = node[1].asInt();
        const int y0 = node[2].asInt();
        const int y1 = node[3].asInt();
        const int z0 = node[4].asInt();
        const int z1 = node[5].asInt();
        const IPPVector3DInt lower = {{x0, y0, z0}};
        const IPPVector3DInt upper = {{x1, y1, z1}};
        return IPPBox3DInt(lower, upper);
    }
}

static void readBoxes(const Json::Value &node, const size_t dims, std::vector<IPPBox3DInt> &boxes)
{
    for (size_t i = 0; i < node.size(); ++i) {
        const Json::Value &entry = node[(int) i];
        boxes.emplace_back(readBox(entry, dims));
    }
}

static void parseCompositions(const Json::Value &compositions,
                              IPPConfig &config,
                              std::map<std::string, Composition *> &compNameToComp)
{
    config.compositions.resize(compositions.size());

    for (size_t iComp = 0; iComp < compositions.size(); ++iComp) {
        Composition &composition = config.compositions[iComp];
        const Json::Value &comp = compositions[(int) iComp];

        composition.name = comp["name"].asString();
        composition.temp = comp.get("temp", 25.0).asDouble();
        composition.pressure = comp.get("pressure", 1.0).asDouble();
        composition.pH = Composition::ScalarAndString(comp.get("pH", 7.0).asDouble(),
                                                      comp.get("pH_extension", "charge").asString());

        composition.pe = Composition::ScalarAndString(comp.get("pe", 4.0).asDouble(),
                                                      comp.get("pe_extension", "").asString());

        composition.units = comp.get("units", "mol/L").asString();
        composition.inertVolFraction = comp.get("inert_fraction", 0.0).asDouble();

        if (composition.inertVolFraction < 0.0 || composition.inertVolFraction > 1.0) {
            throw std::runtime_error("inert volume fraction for component " + composition.name
                                     + " out of range: "
                                     + std::to_string(composition.inertVolFraction));
        }

        compNameToComp[composition.name] = &composition;

        // reading solution concentration
        const Json::Value &element_concentrations = comp["element_concentrations"];
        composition.elemConcVec.resize(element_concentrations.size());

        for (size_t iElem = 0; iElem < element_concentrations.size(); ++iElem) {
            Composition::ElementConcentration &elemConc = composition.elemConcVec[iElem];

            const Json::Value &elem = element_concentrations[(int) iElem];
            elemConc.elemSpecies.elementName = elem["name"].asString();
            elemConc.elemSpecies.speciesName = elem["species"].asString();
            elemConc.conc = elem["conc"].asDouble();
        }

        // reading phases
        const Json::Value &phases = comp["phases"];
        composition.phases.resize(phases.size());

        for (size_t iPhase = 0; iPhase < phases.size(); ++iPhase) {
            Composition::PhaseDefinition &phaseAmount = composition.phases[iPhase];

            const Json::Value &phase = phases[(int) iPhase];
            phaseAmount.name = phase["name"].asString();
            phaseAmount.amount = phase["amount"].asDouble();
        }
    }
}

static IPPConfig::TransportOnlyFlag parseTransportOnlyFlag(const Json::Value &root)
{
    const std::string flag = root.get("transport_only_flag", "").asString();

    if (flag == "init_chem_only") {
        return IPPConfig::TOF_InitChemOnly;
    } else if (flag == "fake_chem") {
        return IPPConfig::TOF_FakeChem;
    } else if (flag.empty() == false) {
        throw std::runtime_error("unknown transport only flag: " + flag);
    } else {
        return IPPConfig::TOF_FullChemistry;
    }
}

static IPPConfigPtr read(int argc,
                         char *argv[],
                         const boost::program_options::variables_map &vm,
                         const std::string &fileName,
                         BCGenerator &bcGen,
                         SaturationIndexCalcFactory &si_calcFactory,
                         PorosityModelFactory &porosModelFactory,
                         MultiScaleDiffusionFactory &multiScaleFactory)
{
    IPPConfigPtr configPtr = std::make_shared<IPPConfig>(argc, argv);

    IPPConfig &config = *configPtr;

    if (vm.count("verbose")) {
        config.verbose = true;
    } else {
        config.verbose = false;
    }

    // pcout << "reading: " << fileName << std::endl;

    const boost::filesystem::path filePath(fileName);
    const boost::filesystem::path fn = filePath.filename();
    const boost::filesystem::path base = fn.stem();
    const std::string fallbackOutput = base.string();

    std::ifstream file(fileName);
    if (file.is_open() == false) {
        throw std::runtime_error("Could not open file: " + fileName);
    }

    Json::Value root;

    {
        const ScopedDisableFloatingPointException fpe;
        file >> root;
    }

    config.outputDir = root.get("outputDir", fallbackOutput).asString();
    config.maxIter = root.get("maxIter", 0).asUInt64();
    config.maxTime = root.get("maxTime", 0).asDouble();
    config.outputFrequency = root.get("outputFrequency", 1).asUInt64();
    config.outletFluxFrequency = root.get("outletFluxFrequency", 0).asUInt64();
    config.checkpointingTime = root.get("checkpointingTime", std::numeric_limits<double>::max())
                                   .asDouble();
    config.checkpointingFreq = root.get("checkpointingFreq", std::numeric_limits<size_t>::max())
                                   .asUInt64();

    config.latticeBoltzmannCollType = getCollisionType(root);
    config.diffusionCoefReference = root.get("diff_coef_reference", 1.0E-9).asDouble();
    config.tauRef = root.get("tau_reference", 1.0).asDouble();
    config.porosRef = root.get("porosity_reference", 1.0).asDouble();
    config.porosLow = root.get("porosity_low", 1.0).asDouble();

    config.isChemistryEnabled = root.get("enable_chemistry", true).asBool();
    config.transportOnlyFlag = parseTransportOnlyFlag(root);
    config.chemistryTolerance = root.get("chemistry_tolerance", 0.0).asDouble();

#ifdef DEBUG_OUTPUT
    config.write_chem = root.get("write_chem", true).asBool();
#else
    config.write_chem = root.get("write_chem", false).asBool();
#endif

    const std::string defaultDBPath = "databases";
    config.dataBasePath = root.get("database_path", defaultDBPath).asString();

    const std::string defaultDB = "phreeqc.dat";
    config.dataBaseName = root.get("database", defaultDB).asString();

    if (root.isMember("database_extension")) {
        const Json::Value &extNode = root["database_extension"];
        const std::string dbFileName = config.dataBasePath + "/" + extNode.asString();

        std::ifstream dbFile(dbFileName);

        if (dbFile.is_open() == false) {
            throw std::runtime_error("Could not open file: " + dbFileName);
        }

        Json::Value dbRoot;
        {
            ScopedDisableFloatingPointException fpe;
            dbFile >> dbRoot;
        }

        IPPConfig::DataBaseExtension &ext = config.dataBaseExt;

        const Json::Value &phases = dbRoot["phases"];
        for (size_t iPhase = 0; iPhase < phases.size(); ++iPhase) {
            const Json::Value &phase = phases[(int) iPhase];
            const std::string phaseName = phase["name"].asString();

            ext.phaseNameToInfos[phaseName] = PhaseInfo();
            PhaseInfo &phaseInfo = ext.phaseNameToInfos[phaseName];

            phaseInfo.molarVolume = phase.get("molarVolume", -1.0).asDouble();
        }
    }

    if (root.isMember("permeable_phases")) {
        const Json::Value &permNode = root["permeable_phases"];

        IPPConfig::DataBaseExtension &ext = config.dataBaseExt;
        for (size_t iPhase = 0; iPhase < permNode.size(); ++iPhase) {
            const Json::Value &phase = permNode[(int) iPhase];
            const std::string phaseName = phase.asString();
            PhaseInfo &phaseInfo = ext.phaseNameToInfos[phaseName];
            phaseInfo.isPermeable = true;
        }
    }

    if (root.isMember("chemistry_prediction")) {
        const Json::Value &node = root["chemistry_prediction"];
        IPPCheck::assertCheck(node.isMember("lookup_file"));
        config.chemPredict.filename = node["lookup_file"].asString();
        config.chemPredict.forceRecalcFreq = node.get("force_recalc_freq", 0).asUInt();
        config.chemPredict.enabled = true;
    }

    // reading solid solutions
    const Json::Value &solidSolutions = root["solid_solutions"];

    for (size_t iSS = 0; iSS < solidSolutions.size(); ++iSS) {
        const Json::Value &solidSolution = solidSolutions[(int) iSS];

        config.solidSolutionVec.push_back(IPPConfig::SolidSolutionDefinition());
        IPPConfig::SolidSolutionDefinition &ssDef = config.solidSolutionVec.back();
        ssDef.name = solidSolution["name"].asString();

        std::vector<std::string> &comps = ssDef.phaseNames;
        const Json::Value &ssComponents = solidSolution["components"];
        for (size_t iComp = 0; iComp < ssComponents.size(); ++iComp) {
            const Json::Value &comp = ssComponents[(int) iComp];
            const std::string &compName = comp.asString();
            comps.push_back(compName);
        }
    }

    // reading composition informations
    std::map<std::string, Composition *> compNameToComp;

    config.isDegradedDiffCalc = false;

    if (root.isMember("compositions")) {
        parseCompositions(root["compositions"], config, compNameToComp);
    } else if (root.isMember("diffusion_degraded")) {
        config.isDegradedDiffCalc = true;
    } else {
        throw std::runtime_error("neither compositions nor diffusion_degraded defined");
    }

    bcGen.setCompositionCookie(&compNameToComp);

    // reading lattice dimensions
    {
        const Json::Value &dimNode = root["dimensions"];
        const Json::Value &dimsNode = dimNode["dims"];
        config.nx = dimsNode.get(0u, 1).asUInt();
        config.ny = dimsNode.get(1u, 1).asUInt();
        config.nz = dimsNode.get(2u, 1).asUInt();

        if (config.nz > 1) {
            config.is3DSimulation = true;
        } else {
            config.is3DSimulation = false;
        }

        const Json::Value &offset = dimNode["offset"];
        config.offset = {{offset.get(0u, 0.0).asDouble(),
                          offset.get(1u, 0.0).asDouble(),
                          offset.get(2u, 0.0).asDouble()}};

        config.spatialResolution = dimNode.get("spatial_resolution", 1.0).asDouble();
        config.distTransPorosThresh = dimNode.get("dist_trans_porosity_threshold", 0.5).asDouble();

        const Json::Value &bc = dimNode["boundary_conditions"];
        const ConfigDataWrapper confWrapper(bc);

        const size_t nDims = config.is3DSimulation ? 3 : 2;

        bcGen.setDim(nDims);
        bcGen.setCookie(&confWrapper);
        config.boundaryConditions.reset(bcGen.generate());

        readBoxes(dimNode["chemistry_disabled"], nDims, config.chemDisabledBoxes);

        const Json::Value enabledChemNode = dimNode["trans_sync_enabled"];
        if (enabledChemNode.isArray()) {
            config.transportSyncEnabledBox = readBox(enabledChemNode, nDims);
        } else {
            // whe not defined update all cells in transport module
            const int min = std::numeric_limits<int>::lowest();
            const int max = std::numeric_limits<int>::max();
            config.transportSyncEnabledBox = IPPBox3DInt({min, min, min}, {max, max, max});
        }

        readBoxes(dimNode["decomp"], nDims, config.decomp);
    }
    const int nxyz = config.nx * config.ny * config.nz;

    if (root.isObject() && root.isMember("advection")) {
        const Json::Value &advNode = root["advection"];
        const Json::Value &flowVecNode = advNode["flow_vector"];
        const IPPVector3D flowVec = {{flowVecNode.get(0u, 0.0).asDouble(),
                                      flowVecNode.get(1u, 0.0).asDouble(),
                                      flowVecNode.get(2u, 0.0).asDouble()}};

        const double fluidDensity = advNode.get("fluid_density", 1.0).asDouble();
        config.flowFunc = AbstractFlowFunctorPtr(
            new InitFunctors::ConstantVelocityAndDensity(flowVec, fluidDensity));
    }

    IDtoComposition idToComposition;

    if (config.isDegradedDiffCalc) {
        std::string phases;
        if (vm.count("phases")) {
            phases = vm["phases"].as<std::string>();
        }

        std::string comp;
        if (vm.count("molarity")) {
            comp = vm["molarity"].as<std::string>();
        }

        if (phases.empty() || comp.empty()) {
            // there must be at least the phases and comp VTK output files
            throw std::runtime_error(
                "degraded VTK input set must defined as parameters: <phases.vti> <comp.vti>");
        }

        const Json::Value degNode = root["diffusion_degraded"];

        const std::string inertFracMap = degNode.get("inert_frac", "").asString();

        const Json::Value &offset = degNode["origin"];
        IPPVector3DInt degOffset = {
            {offset.get(0u, 0).asInt(), offset.get(1u, 0).asInt(), offset.get(2u, 0).asInt()}};

        Composition compDefault;
        compDefault.units = "mol/L";
        compDefault.pH.value = degNode.get("pH", 7.0).asDouble();
        compDefault.temp = degNode.get("temp", 25).asDouble();

        VTKImageReader::readDegrVTKs(phases, comp, inertFracMap, compDefault, degOffset, config);
    } else {
        // reading ID mapping
        const Json::Value &idMapping = root["ID_mapping"];
        for (size_t iMapping = 0; iMapping < idMapping.size(); ++iMapping) {
            const Json::Value &value = idMapping[(int) iMapping];
            const Json::Value::Members keys = value.getMemberNames();
            const std::string &compositionName = keys.front();
            const size_t ID = value[compositionName].asInt();

            IDtoComposition::const_iterator it = idToComposition.find(ID);
            if (it != idToComposition.end()) {
                const std::string errStr = "ID mapping ambigous. ID -> composition name: "
                                           + std::to_string(ID) + " -> " + compositionName + "vs."
                                           + std::to_string(it->first) + " -> " + it->second->name;
                throw std::runtime_error(errStr);
            }

            if (compNameToComp.find(compositionName) == compNameToComp.end()) {
                throw std::runtime_error("Composition not found in input file: "
                                         + std::to_string(ID) + " -> " + compositionName);
            }

            Composition *compPtr = compNameToComp[compositionName];
            idToComposition[ID] = compPtr;
        }

        const Json::Value geometryNode = root["geometry"];
        readGeometry(geometryNode, idToComposition, config);
    }

    // check input data consistency
    for (const Domain &domain : config.domains) {
        for (int cellID : domain.cells) {
            if (cellID < 0 || cellID >= nxyz) {
                throw std::runtime_error(
                    "geometry definition is larger than lattice dimensions provide");
            }
        }
    }

    const Json::Value &whiteList = root["phase_white_list"];
    for (size_t i = 0; i < whiteList.size(); ++i) {
        const Json::Value &entry = whiteList[(int) i];
        const std::string phaseName = entry.asString();
        config.whiteListedPhases.insert(phaseName);
    }

    const Json::Value &blackList = root["phase_black_list"];
    for (size_t i = 0; i < blackList.size(); ++i) {
        const Json::Value &entry = blackList[(int) i];
        const std::string phaseName = entry.asString();
        config.blackListedPhases.insert(phaseName);
    }

    const Json::Value &waterLimitedList = root["water_limited_phases"];
    for (size_t i = 0; i < waterLimitedList.size(); ++i) {
        const Json::Value &entry = waterLimitedList[(int) i];
        const std::string phaseName = entry.asString();
        config.waterLimitedPhases.insert(phaseName);
        config.kineticPhases.insert(phaseName);
    }

    const Json::Value &kineticPhasesNode = root["kinetic_phases"];
    for (Json::ArrayIndex i = 0; i < kineticPhasesNode.size(); ++i) {
        const Json::Value phaseNode = kineticPhasesNode[i];
        const std::string phaseName = phaseNode.asString();
        config.kineticPhases.insert(phaseName);
    }

    const Json::Value &phaseConf = root["phase_diss_prec_config"];

    PhaseNameToDissPrecipBehaviour *dissPrec = new PhaseNameToDissPrecipBehaviour;
    config.dissPrecBehav.reset(dissPrec);

    for (size_t i = 0; i < phaseConf.size(); ++i) {
        const Json::Value &entry = phaseConf[(int) i];
        const std::string phaseName = entry["name"].asString();
        const std::string value = entry["value"].asString();

        DissPrecipBehaviour data = DPB_Normal;
        if (value == "precipitate_only") {
            data = DPB_PrecOnly;
        } else if (value == "dissolve_only") {
            data = DPB_DissOnly;
        } else {
            throw std::runtime_error("Unknown value for phase config: " + phaseName + "\t" + value);
        }

        dissPrec->addData(phaseName, data);
    }

    if (root.isMember("porosity_model")) {
        const Json::Value &porosityNode = root["porosity_model"];
        const ConfigDataWrapper confWrapper(porosityNode);
        porosModelFactory.setCookie(&confWrapper);
        config.porosCalc.reset(porosModelFactory.generate());
    }

    if (root.isMember("multiscale_func")) {
        const Json::Value &multiscaleNode = root["multiscale_func"];
        const ConfigDataWrapper confWrapper(multiscaleNode);
        multiScaleFactory.setCookie(&confWrapper);
        config.multiScaleDiffusionCalc.reset(multiScaleFactory.generate());
    }

    if (root.isMember("saturation_index_func")) {
        si_calcFactory.setSpatialResolution(config.spatialResolution);
        const Json::Value &siNode = root["saturation_index_func"];
        const ConfigDataWrapper confWrapper(siNode);
        si_calcFactory.setCookie(&confWrapper);
        si_calcFactory.setDissPrecipitationBehaviour(dissPrec);
        config.saturationIndexCalc.reset(si_calcFactory.generate());
    }

    const double porosThreshDefault = 0.001;
    AbstractDissolveOnlyCalc *calc = nullptr;
    if (root.isMember("dissolve_only_func")) {
        const Json::Value &doNode = root["dissolve_only_func"];
        calc = parseDissolveOnlyFunc(config.is3DSimulation, doNode, porosThreshDefault);
    } else {
        std::cout << "WARNING: No dissolve only function defined, fallback to simple minimum "
                     "porosity technique. Thresh value: "
                  << porosThreshDefault << std::endl;
        calc = new BelowThreshDissolveOnly(porosThreshDefault);
    }

    config.dissolveOnlyFunc.reset(calc);

    parseResults(root, config);

    config.init();

    return configPtr;
}

static ConfigReader::Factories checkFactories(const ConfigReader::Factories &factories)
{
    ConfigReader::Factories checked = factories;

    if (checked.bcGen == nullptr) {
        checked.bcGen.reset(new NonReturningBCGen);
    }

    if (checked.si_calcFactory == nullptr) {
        checked.si_calcFactory.reset(new NonReturningSICalcFactory);
    }

    if (checked.porosModelFactory == nullptr) {
        checked.porosModelFactory.reset(new NonReturningPorosityModelFactory);
    }

    if (checked.multiScaleFactory == nullptr) {
        checked.multiScaleFactory.reset(new NonReturningMultiScaleDiffusionFactory);
    }

    return checked;
}

IPPConfigPtr ConfigReader::read(int argc,
                                char *argv[],
                                const boost::program_options::variables_map &vm,
                                const std::string &fileName,
                                const Factories &factories)
{
    const Factories checked = checkFactories(factories);

    return read(argc,
                argv,
                vm,
                fileName,
                *checked.bcGen,
                *checked.si_calcFactory,
                *checked.porosModelFactory,
                *checked.multiScaleFactory);
}

IPPConfigPtr ConfigReader::read(int argc,
                                char *argv[],
                                const boost::program_options::variables_map &vm,
                                const std::string &fileName)
{
    const Factories factories;
    return read(argc, argv, vm, fileName, factories);
}

} // namespace IPP
