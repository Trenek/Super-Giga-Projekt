# ── project ──────────────────────────────────────────────────────────────────
PROGS   = output/cubicIkeda
OTHERS  = computations/integrateFuncs computations/toOdeFuncs

# ── capd ─────────────────────────────────────────────────────────────────────
CAPDBINDIR =
CAPDCXX    := $(shell $(CAPDBINDIR)capd-config --variable=capd_cxx)
CAPDFLAGS   = $(shell $(CAPDBINDIR)capd-config --cflags)
CAPDLIBS    = $(shell $(CAPDBINDIR)capd-config --libs)

# ── compiler ─────────────────────────────────────────────────────────────────
CXXFLAGS += $(CAPDFLAGS) -O2 -std=c++17 -Icomputations

# ── internals ────────────────────────────────────────────────────────────────
OBJDIR     = .obj/
MAIN_OBJ   = $(OBJDIR)main.o
OTHERS_OBJ = $(OTHERS:%=$(OBJDIR)%.o)
OBJ_FILES  = $(MAIN_OBJ) $(OTHERS_OBJ)

# ── targets ──────────────────────────────────────────────────────────────────
.PHONY: all run plot clean

all: $(PROGS)

$(PROGS): $(MAIN_OBJ) $(OTHERS_OBJ) | output
	$(CAPDCXX) -o $@ $(MAIN_OBJ) $(OTHERS_OBJ) $(CAPDLIBS)

run: all
	./output/cubicIkeda

plot: run
	venv/bin/python3 plot_traj.py

output:
	mkdir -p output

# ── compile + dependency generation ──────────────────────────────────────────
-include $(OBJ_FILES:%=%.d)

$(OBJ_FILES): $(OBJDIR)%.o : %.cpp
	@mkdir -p $(dir $@)
	$(CAPDCXX) $(CXXFLAGS) -MT $@ -MD -MP -MF $(@:%=%.d) -c -o $@ $<

# ── clean ────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OBJDIR) $(PROGS)