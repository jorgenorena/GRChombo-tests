#ifndef RADIATIONACOUSTICLEVEL_HPP_
#define RADIATIONACOUSTICLEVEL_HPP_
#include "AcousticMatter.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class RadiationAcousticLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<RadiationAcousticLevel>;
    using GRAMRLevel::GRAMRLevel;
    WallRadiation matter() const
    {
        return WallRadiation(Potential(m_p.potential_params), m_p.radiation_kappa);
    }
    void refresh_diagnostics();
    void specificAdvance() override;
    void initialData() override;
    void postRestart() override { refresh_diagnostics(); }
    void specificEvalRHS(GRLevelData &soln, GRLevelData &rhs, double time) override;
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
