enum ClusteringCriterion {
  WARD_D = 1,
  SINGLE = 2,
  COMPLETE = 3,
  AVERAGE = 4,
  MCQUITTY = 5,
  MEDIAN = 6,
  CENTROID = 7,
  WARD_D2 = 8
};

/* A signed integer type to store the index of a vector
 * I would use R_xlen_t, but I would not be able to return an object with those
 * indices (I'd have to coerce them to... doubles?)
 * It's a bit confusing but I want to decide later what to use, so I'll just
 * typedef this, use an "int" and say "long vectors not supported".
 */
typedef int idx_t;

int hclust(idx_t n, int iopt, idx_t *ia, idx_t *ib, double *crit, double *membr,
           idx_t *nn, double *disnn, double *diss, int *flag);

int hcass2(idx_t n, const idx_t *ia, const idx_t *ib, idx_t *iorder, idx_t *iia, idx_t *iib);

