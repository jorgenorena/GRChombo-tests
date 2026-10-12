GRCHOMBO_SOURCE = ../../Source
ebase := RadiationMatterTest
LibNames := BoxTools
cxxcppflags += -I$(GRCHOMBO_SOURCE)/GRChomboCore
src_dirs := $(GRCHOMBO_SOURCE)/utils $(GRCHOMBO_SOURCE)/simd \
            $(GRCHOMBO_SOURCE)/CCZ4 $(GRCHOMBO_SOURCE)/Matter \
            $(GRCHOMBO_SOURCE)/BoxUtils
include $(CHOMBO_HOME)/mk/Make.test
$(objects): $(wildcard *.hpp) $(wildcard $(GRCHOMBO_SOURCE)/Matter/*.hpp) \
            ../../Examples/RadiationAcoustic/UserVariables.hpp \
            ../../Examples/RadiationAcoustic/DiagnosticVariables.hpp
