### Microbial turnover across Arctic and North Atlantic environmental gradients

Marine microorganisms sustain major ocean biogeochemical cycles, yet their
community composition varies across water masses, depth layers and
physicochemical gradients [1,2]. We reanalysed ribosomal-RNA miTAG community
profiles and matched environmental metadata from the Tara Oceans expedition,
which sampled planktonic communities across the global ocean between 2009 and
2013 [2,3]. The source dataset contained 180 prokaryote-enriched metagenomic
samples. We selected 38 Arctic Ocean samples to represent an environmentally
extreme and comparatively understudied polar ecosystem relative to
lower-latitude oceans. We selected 24 North Atlantic Ocean samples because
their spatial distribution broadly traces major circulation pathways, which
can delimit plankton genomic provinces [4]. Together, the selected samples span
surface, deep-chlorophyll-maximum, mixed and mesopelagic waters. The mapped
currents provide geographical context only; current transport was not included
as a covariate in the present analysis.

We treated each water sample as a local community and each ocean region as its
regional pool. HRIC alpha and gamma diversity were calculated at the sample and
region levels, respectively, and sample-level turnover was defined as gamma
minus alpha. Mean alpha diversity was 0.270 in the Arctic and 0.368 in the
North Atlantic, whereas regional gamma diversity was 0.379 and 0.503,
respectively. The corresponding additive beta diversities, equal to mean
turnover within each region, were 0.109 and 0.135. These values describe the
two selected regional pools; they are not a formal test of an Arctic--North
Atlantic difference. Ordering samples from the smallest to the largest
turnover linked each local diversity value to the taxonomic composition shown
directly below it.

Two environmental trends were consistent across regions. Turnover declined
with depth in both the Arctic (Spearman's rho = -0.68, P = 5.1 x 10^-6, n =
36) and North Atlantic (rho = -0.60, P = 0.0017, n = 24), and increased with
chlorophyll a (Arctic, rho = 0.64, P = 1.7 x 10^-5, n = 38; North Atlantic,
rho = 0.68, P = 2.5 x 10^-4, n = 24). The Arctic gradient also included higher
turnover at greater oxygen concentration (rho = 0.61, P = 4.6 x 10^-5, n = 38)
and lower turnover at greater phosphate (rho = -0.59, P = 0.0010, n = 28) and
nitrate concentrations (rho = -0.66, P = 1.2 x 10^-4, n = 28). Arctic
turnover was not associated with temperature (rho = -0.08, P = 0.65, n = 38).
In the North Atlantic, turnover increased with temperature (rho = 0.44, P =
0.033, n = 24) and declined with nitrate (rho = -0.48, P = 0.023, n = 22).
Associations with oxygen (P = 0.051) and phosphate (P = 0.33) were not
supported at the 0.05 level. These are exploratory, univariable correlations
with unadjusted two-sided P values. The environmental variables covary in the
water column, so their coefficients should be interpreted as gradients rather
than independent effects.

The stacked profiles showed that samples with similar turnover could contain
different dominant taxa, and that the taxonomic signature of the environmental
gradient differed between regions. We used nitrate for the taxon-resolved
panel because, among the six displayed covariates, it yielded the largest
number of significant univariable associations with HRIC-transformed taxon
coordinates. This choice was data-driven, so the taxon-level results are
exploratory. Among the 19 most abundant taxa in each region, 14 Arctic and 13
North Atlantic taxa were associated with nitrate after Benjamini--Hochberg
adjustment within region (q < 0.05). In the Arctic, the HRIC coordinates of
SAR86, SAR11 and Marinimicrobia increased with nitrate, whereas those of
Formosa, Cryomorphaceae, NS5, Polaribacter, SAR92, Marinoscillum and
Flavobacteriaceae decreased. In the North Atlantic, SAR324, SAR202 and
Marinimicrobia increased with nitrate, whereas NS4, OM1, Gammaproteobacteria,
Prochlorococcus, NS5, SAR86 and Marinoscillum decreased. Marinimicrobia showed
the same positive association in both regions, and NS5, chloroplast-assigned
sequences and Marinoscillum showed negative associations in both. By contrast,
SAR86 changed sign between regions. These coefficients describe changes in
HRIC taxon coordinates, not changes in relative abundance alone, and do not by
themselves identify causal drivers of turnover.

Together, the regional analyses identify a shared depth--productivity axis
along which microbial turnover changed in both oceans, superimposed on
region-specific temperature, oxygen and nutrient patterns. The taxon-resolved
results further show that a similar nitrate--turnover trend need not arise from
the same community reorganization in each basin. Because the samples are
observational, unevenly distributed among water-column layers and exposed to
correlated environmental gradients, these results support environmental
filtering as a hypothesis rather than a causal conclusion. They nevertheless
demonstrate how local diversity, regional diversity, turnover and individual
taxon coordinates can be followed in one coherent analysis.

**Figure caption. Environmental gradients accompany microbial turnover across
Arctic and North Atlantic waters.** **a--d,** Arctic Ocean communities (n =
38). **e--h,** North Atlantic Ocean communities (n = 24). **a,e,** Sampling
locations, shown as black crosses. Schematic warm and cold currents in **e**
provide geographical context and were not included in the statistical models.
**b,f,** Sample-level alpha diversity calculated using `HRIC::SHalpha`.
Samples are ordered from the smallest to the largest turnover, defined as gamma
minus alpha; dashed red lines show region-specific gamma diversity calculated
using `HRIC::SHgamma`. **c,g,** Relative abundance of the 19 taxa with the
largest mean relative abundance in each region, with all remaining taxa grouped
as Others. Bars follow the sample order in **b** and **f**. Asterisks identify
taxa whose HRIC-transformed coordinate was associated with nitrate in
univariable linear models at Benjamini--Hochberg-adjusted q < 0.05, with
adjustment across the 19 tested taxa within each region. **d,h,** Turnover in
relation to depth, temperature, oxygen, chlorophyll a, phosphate and nitrate.
Red lines and shaded areas are ordinary least-squares visual summaries with 95%
confidence intervals. Annotations report Spearman's rho, unadjusted two-sided P
values and complete-case sample sizes.

#### References

1. Sunagawa, S. et al. Structure and function of the global ocean microbiome.
   *Science* **348**, 1261359 (2015).
2. Salazar, G. et al. Gene expression changes and community turnover
   differentially shape the global ocean metatranscriptome. *Cell* **179**,
   1068--1083.e21 (2019).
3. Pesant, S. et al. Open science resources for the discovery and analysis of
   Tara Oceans data. *Scientific Data* **2**, 150023 (2015).
4. Richter, D. J. et al. Genomic evidence for global ocean plankton
   biogeography shaped by large-scale current systems. *eLife* **11**, e78129
   (2022).
