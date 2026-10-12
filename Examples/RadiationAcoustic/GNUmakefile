GRCHOMBO_SOURCE = ../../Source
ebase := Main_RadiationAcoustic
LibNames := AMRTimeDependent AMRTools BoxTools
src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd \
            $(GRCHOMBO_SOURCE)/CCZ4 \
            $(GRCHOMBO_SOURCE)/Matter \
            $(GRCHOMBO_SOURCE)/BoxUtils \
            $(GRCHOMBO_SOURCE)/GRChomboCore \
            $(GRCHOMBO_SOURCE)/AMRInterpolator \
            $(GRCHOMBO_SOURCE)/ApparentHorizonFinder
include $(CHOMBO_HOME)/mk/Make.test
$(objects): $(wildcard *.hpp) $(wildcard ../DW_flat/*Potential.hpp) \
            ../RadiationFLRW/FLRWDiagnostics.hpp \
            $(wildcard $(GRCHOMBO_SOURCE)/Matter/*.hpp) \
            $(wildcard $(GRCHOMBO_SOURCE)/CCZ4/*.hpp)
