#ifndef RADIATIONHANDOFFLEVEL_HPP_
#define RADIATIONHANDOFFLEVEL_HPP_
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "../RadiationFLRW/FLRWDiagnostics.hpp"

class RadiationHandoffLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<RadiationHandoffLevel>;
    using GRAMRLevel::GRAMRLevel;
    HDF5HeaderData m_initial_header;
    WallRadiation matter() const
    { return WallRadiation(Potential(m_p.potential_params), m_p.background.kappa, true); }
    void refresh_diagnostics();
    void initialData() override { MayDay::Error("RadiationHandoff requires a GRTresna restart file"); }
    void postRestart() override;
    void specificEvalRHS(GRLevelData &, GRLevelData &, double) override
    { MayDay::Error("Pass D never advances the imported initial slice"); }
    void computeTaggingCriterion(FArrayBox &criterion, const FArrayBox &, const FArrayBox &) override
    { criterion.setVal(0.); }
    void prePlotLevel() override { refresh_diagnostics(); }
    void writeCheckpointHeader(HDF5Handle &handle) const override;
    void readCheckpointHeader(HDF5Handle &handle) override;
  public:
    void specificPostTimeStep() override;
};
#endif
