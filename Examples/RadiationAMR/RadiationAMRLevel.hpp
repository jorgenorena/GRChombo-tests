#ifndef RADIATIONAMRLEVEL_HPP_
#define RADIATIONAMRLEVEL_HPP_
#include "AcousticMatter.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include <array>

class RadiationAMRLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<RadiationAMRLevel>;
    using GRAMRLevel::GRAMRLevel;
    WallRadiation matter() const
    { return WallRadiation(Potential(m_p.potential_params), m_p.radiation_kappa); }
    std::array<double, 5> m_before_transfer{};
    void check_radiation(GRLevelData &state);
    void refresh_diagnostics();
    void refresh_hierarchy();
    // Q, energy, minimum timelike margin, active levels, fine-cell count.
    std::array<double, 5> transfer_state();
    void write_layout();
    void initialData() override;
    void specificAdvance() override;
    void specificEvalRHS(GRLevelData &soln, GRLevelData &rhs, double time) override;
    void preTagCells() override { refresh_diagnostics(); }
    void computeTaggingCriterion(FArrayBox &criterion, const FArrayBox &state,
                                 const FArrayBox &diagnostics) override;
    void preRegrid(int base, const Vector<Vector<Box>> &new_grids) override;
    void specificPostRegrid(int base) override;
#ifdef CH_USE_HDF5
    void postRestart() override { refresh_diagnostics(); }
    void preCheckpointLevel() override { refresh_diagnostics(); }
    void prePlotLevel() override { refresh_diagnostics(); }
    void writeCheckpointHeader(HDF5Handle &handle) const override;
    void readCheckpointHeader(HDF5Handle &handle) override;
#endif
  public:
    void specificPostTimeStep() override;
};
#endif
