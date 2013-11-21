EXEDIR := scripts
OBJDIR := bin
SRCDIR := src
INCDIR := inc
MAKEDIR := bin

CXX := $(shell root-config --cxx)
CXXFLAGS := -isystem $(shell root-config --incdir) -Wall -Wextra -pedantic -Wshadow $(shell root-config --cflags) -O2 -I $(INCDIR)
LD := $(shell root-config --ld)
LDFLAGS := $(shell root-config --ldflags)
LDLIBS := $(shell root-config --libs) -lMinuit

vpath %.cpp $(SRCDIR)
vpath %.hpp $(INCDIR)
vpath %.o $(OBJDIR)
vpath %.so $(OBJDIR)
vpath %.exe $(EXEDIR)
vpath %.d $(MAKEDIR)

# Add new executables to this list
all: make_plots.exe skim_file.exe stack_histos.exe draw_abcd_ratio_plots.exe make_sig_plots.exe calc_abcd.exe count_specific_mass_events.exe draw_npv_plot.exe make_cutflow_table.exe calc_abcd_new.exe piechart.exe

# List any object files your executable need to be linked with
$(EXEDIR)/draw_npv_plot.exe: draw_npv_plot.o pu_constants.o
$(EXEDIR)/count_specific_mass_events.exe: count_specific_mass_events.o
$(EXEDIR)/calc_abcd.exe: calc_abcd.o weights.o
$(EXEDIR)/calc_abcd_new.exe: calc_abcd_new.o weights.o math.o abcd_calculator.o abcd_count.o
$(EXEDIR)/piechart.exe: piechart.o
$(EXEDIR)/make_sig_plots.exe: make_sig_plots.o
$(EXEDIR)/generate_cfa_class.exe: generate_cfa_class.o
$(EXEDIR)/stack_histos.exe: stack_histos.o
$(EXEDIR)/draw_abcd_ratio_plots.exe: draw_abcd_ratio_plots.o
$(EXEDIR)/skim_file.exe: skim_file.o lib_jet_met_objects.so event_handler.o event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o
$(EXEDIR)/make_plots.exe: make_plots.o lib_jet_met_objects.so event_handler.o event_number.o b_jet.o math.o pu_constants.o timer.o cfa.o weights.o
$(EXEDIR)/make_cutflow_table.exe: make_cutflow_table.o cutflow.o

-include $(addsuffix .d,$(addprefix $(MAKEDIR)/,$(notdir $(basename $(wildcard $(SRCDIR)/*.cpp)))))
-include $(MAKEDIR)/cfa.d

$(MAKEDIR)/%.d: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM -MG -MF $@ $< 
	sed -i'' 's#$*.o#$(OBJDIR)/$*.o $(MAKEDIR)/$*.d#g' $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# This is a bit ugly. Shouldn't need the dependency explicitly.
$(EXEDIR)/%.exe: $(OBJDIR)/%.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# cfa.cpp and cfa.hpp need special treatment. Probably cleaner ways to do this.
$(SRCDIR)/cfa.cpp $(INCDIR)/cfa.hpp: dummy_cfa.all
.SECONDARY: dummy_cfa.all
dummy_cfa.all: $(EXEDIR)/generate_cfa_class.exe example_root_file.root
	./$< $(word 2,$^)
.PRECIOUS: generate_cfa_class.o

.DELETE_ON_ERROR:

.PHONY: clean

clean:
	-rm -rf $(EXEDIR)/*.exe $(OBJDIR)/*.o $(MAKEDIR)/*.d $(SRCDIR)/cfa.cpp $(INCDIR)/cfa.hpp *.exe *.o *.d
	./scripts/remove_backups.sh
