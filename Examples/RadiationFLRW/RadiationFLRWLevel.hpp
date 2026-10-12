#ifndef RADIATIONFLRWLEVEL_HPP_
#define RADIATIONFLRWLEVEL_HPP_
#include "DefaultLevelFactory.hpp"
#include "FLRWDiagnostics.hpp"
#include "GRAMRLevel.hpp"

class RadiationFLRWLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<RadiationFLRWLevel>;
    using GRAMRLevel::GRAMRLevel;
    WallRadiation matter() const
    {
        return WallRadiation(Potential(m_p.potential_params), m_p.background.kappa,
                             m_p.radiation_enabled);
    }
    void refresh_diagnostics();
    void specificAdvance() override;
    void initialData() override;
    void postRestart() override;
    void specificEvalRHS(GRLevelData &soln, GRLevelData &rhs, double time) override;
    void specificUpdateODE(GRLevelData &soln, const GRLevelData &rhs, Real dt) override;
    void preTagCells() override { refresh_diagnostics(); }
    void computeTaggingCriterion(FArrayBox &criterion, const FArrayBox &state,
                                  const FArrayBox &diagnostics) override;
#ifdef CH_USE_HDF5
    void prePlotLevel() override { refresh_diagnostics(); }
    void writeCheckpointHeader(HDF5Handle &handle) const override;
    void readCheckpointHeader(HDF5Handle &handle) override;
#endif
  public:
    void specificPostTimeStep() override;
};
#endif
