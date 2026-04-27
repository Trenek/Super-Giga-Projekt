# ── project ──────────────────────────────────────────────────────────────────
PROGS   = help #periodicOrbit
OTHERS  = computations/integrateFuncs computations/toOdeFuncs computations/newton

# ── CAPD ──────────────────────────────────────────────────────────────────────
CAPDBINDIR =
CAPDCXX   := $(shell $(CAPDBINDIR)capd-config --variable=capd_cxx)
CAPDFLAGS  = $(shell $(CAPDBINDIR)capd-config --cflags)
CAPDLIBS   = $(shell $(CAPDBINDIR)capd-config --libs)

# ── compiler ──────────────────────────────────────────────────────────────────
CXXFLAGS += $(CAPDFLAGS) -O2 -std=c++17 -Icomputations

# ── directories ──────────────────────────────────────────────────────────────
OBJDIR = .obj/

# ── object files ─────────────────────────────────────────────────────────────
OTHERS_OBJ = $(OTHERS:%=$(OBJDIR)%.o)
OBJ_FILES  = $(OTHERS_OBJ) $(PROGS:%=$(OBJDIR)%.o)

.PHONY: all clean run

# ── build all programs ───────────────────────────────────────────────────────
all: $(PROGS:%=output/%)

# ── linking rule (CAPD standard) ─────────────────────────────────────────────
output/%: $(OBJDIR)%.o $(OTHERS_OBJ)
	@mkdir -p output
	$(CAPDCXX) -o $@ $< $(OTHERS_OBJ) $(CAPDLIBS)

# ── dependency tracking ──────────────────────────────────────────────────────
-include $(OBJ_FILES:%=%.d)

# ── compilation rule ──────────────────────────────────────────────────────────
$(OBJ_FILES): $(OBJDIR)%.o : %.cpp
	@mkdir -p $(dir $@)
	$(CAPDCXX) $(CXXFLAGS) -MT $@ -MD -MP -MF $(@:%=%.d) -c -o $@ $<

# ── run ──────────────────────────────────────────────────────────────────────
run: all
# 	./output/periodicOrbit
	./output/help

# ── plot ──────────────────────────────────────────────────────────────────────
plot: run
	venv/bin/python3 plot_traj.py

# ── clean ────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OBJDIR) output