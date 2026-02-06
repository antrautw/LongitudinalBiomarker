# LongitudinalBiomarker
LongitudinalBiomarkerTrajectories+nonlinearInference_V2
Each plot shows three layers: (1) boxplots at each visit summarizing the distribution within each disease, (2) faint within-subject lines 
connecting repeated measures, and (3) a bold disease mean trajectory with a small “slope” label summarizing first-to-last change. For statistics, 
we separate two goals: the main longitudinal inference asks whether diseases differ in how the biomarker changes over time; we use a 
nonlinear mixed model (GAMM) when there are enough visits to support curvature, and automatically fall back to a simpler linear mixed model 
when only two timepoints are available. Separately, the per-visit labels answer a narrower question—whether diseases differ at that specific visit—
and we choose ANOVA only when normality/variance assumptions look reasonable; otherwise we use the more robust Kruskal–Wallis test. This gives us a reviewer-
friendly approach that is both flexible for nonlinear biology and conservative when the data are sparse.

Requires the following variables in numericMeta:
set.label (e.x. b01.127N, b01.127C, etc.)
subject.id (unique per subject identifier)
FirstVLast (contains First, Last, or NA) OR visit.number
Disease OR Group

Slight editing in 
