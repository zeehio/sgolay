/*
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Parameters:                                               C
C                                                            C
C  N                 the number of points being clustered    C
C  DISS(LEN)         dissimilarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, DISNN  vectors of length N, used to store      C
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the dissimilarity assoc. C
C                    with the latter.                        C
C                    MEMBR must be initialized by R to the   C
C                    default of  rep(1, N)                   C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C  Modifications for R: Ross Ihaka, Dec 1996                 C
C                       Fritz Leisch, Jun 2000               C
C  all vars declared:   Martin Maechler, Apr 2001            C
C                                                            C
c- R Bug PR#4195 fixed "along" qclust.c, given in the report C
C- Testing: --> "hclust" in ../../../../tests/reg-tests-1b.R C
C  "ward.D2" (iOpt = 8): Martin Maechler, Mar 2014           C
C                                                            C
C  FLAG not passed as arg to avoid possible                  C
C     C/Fortran inconsistency, May 2019                      C
C                                                            C
C  Translated from Fortran to C by Sergio Oller 2022         C
C  The translation changes ia, ib, and nn, since Fortran is  C
C  1-based indexed and C is 0-based indexed. C stores        C
C  indices in ia, ib and nn with a 0-based criteria          C
C------------------------------------------------------------C
*/

#include <assert.h>
#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>

#include "hclust.h"

#ifndef DBL_MAX
#define DBL_MAX 1E300
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

/* A signed integer type to store the index of a vector
 * I would use R_xlen_t, but I would not be able to return an object with those
 * indices (I'd have to coerce them to... doubles?)
 * It's a bit confusing but I want to decide later what to use, so I'll just
 * typedef this, use an "int" and say "long vectors not supported".
 */
typedef int idx_t;

static inline idx_t ioffst(idx_t n, idx_t i, idx_t j) {
  /*
     Map row I and column J of upper half diagonal symmetric matrix
     onto vector.
   */
  /* prevent idx_t overflows by casting to 64-bit integer */
  int64_t n8 = n, i8 = i, j8 = j;
  return (idx_t)((j8 + i8 * n8 - ((i8+1) * (i8+2)) / 2));
}

int hclust(idx_t n, int iopt, idx_t *ia, idx_t *ib, double *crit, double *membr,
           idx_t *nn, double *disnn, double *diss, int *flag) {
  idx_t npairs = n * (n - 1) / 2;
  /*
   * int ia[n]
   * int ib[n]
   * int nn[n]
   * double crit[n]
   * double membr[n]
   * double disnn[n]
   * double diss[npairs]
   * int flag[n] // could be bool
   */
   idx_t im, jm;
  idx_t i, j, ind, i2, j2, k, ind1, ind2, jj;
  double dmin;
  /* Initializations */
  for (idx_t i = 0; i < n; ++i) {
    // We do not initialize membr in order to be able to restart the algorithm
    // from a cut
    // membr[i] = 1.0;
    flag[i] = 1; /* true */
  }
  if (iopt == WARD_D2) {
    for (idx_t i = 0; i < npairs; ++i) {
      diss[i] = diss[i] * diss[i];
    }
  }

  jj = 0;

  /*
   C
   C  Carry out an agglomeration - first create list of NNs
   C  Note NN and DISNN are the nearest neighbour and its distance
   C  TO THE RIGHT of I.
   *
   */
  for (idx_t i = 0; i < n - 1; ++i) {
    dmin = DBL_MAX;
    idx_t jm = 0;
    for (idx_t j = i + 1; j < n; ++j) {
      ind = ioffst(n, i, j);
      if (dmin > diss[ind]) {
        dmin = diss[ind];
        jm = j;
      }
    }
    // FIXME: if here dmin == DBL_MAX, error, or something
    nn[i] = jm;
    disnn[i] = dmin;
  }

  im = 0;
  jm = 0;
  for (idx_t agl = 0; agl < n - 1; ++agl) {
    // Next, determine least diss. using list of NNs
    dmin = DBL_MAX;
    for (i = 0; i < n - 1; i++) {
      if (flag[i] && disnn[i] < dmin) {
        dmin = disnn[i];
        im = i;
        jm = nn[i];
      }
    }
    // FIXME: if here dmin == DBL_MAX, error or something

    // This allows an agglomeration to be carried out.
    i2 = min(im, jm);
    j2 = max(im, jm);
    ia[agl] = i2;
    ib[agl] = j2;
    if (iopt == WARD_D2)
      dmin = sqrt(dmin);
    crit[agl] = dmin;
    flag[j2] = 0;

    // Update dissimilarities from new cluster.
    dmin = DBL_MAX;
    for (k = 0; k < n; k++) {
      if (flag[k] && k != i2) {
        if (i2 < k) {
          ind1 = ioffst(n, i2, k);
        } else {
          ind1 = ioffst(n, k, i2);
        }
        if (j2 < k) {
          ind2 = ioffst(n, j2, k);
        } else {
          ind2 = ioffst(n, k, j2);
        }
        double d12 = diss[ioffst(n, i2, j2)];

        if (iopt == WARD_D || iopt == WARD_D2) { // WARD'S "D1", or "D2" MINIMUM
                                                 // VARIANCE METHOD - IOPT=8.
          diss[ind1] = (membr[i2] + membr[k]) * diss[ind1] +
                       (membr[j2] + membr[k]) * diss[ind2] - membr[k] * d12;
          diss[ind1] = diss[ind1] / (membr[i2] + membr[j2] + membr[k]);
        } else if (iopt == SINGLE) { // SINGLE LINK METHOD
          diss[ind1] = min(diss[ind1], diss[ind2]);
        } else if (iopt == COMPLETE) {
          diss[ind1] = max(diss[ind1], diss[ind2]);
        } else if (iopt == AVERAGE) { // AVERAGE LINK (OR GROUP AVERAGE) METHOD
          diss[ind1] = (membr[i2] * diss[ind1] + membr[j2] * diss[ind2]) /
                       (membr[i2] + membr[j2]);
        } else if (iopt == MCQUITTY) {
          diss[ind1] = (diss[ind1] + diss[ind2]) / 2;
        } else if (iopt ==
                   MEDIAN) { // MEDIAN (GOWER'S) METHOD aka "Weighted Centroid"
          diss[ind1] = ((diss[ind1] + diss[ind2]) - d12 / 2) / 2;
        } else if (iopt == CENTROID) { // Unweighted CENTROID METHOD
          diss[ind1] = (membr[i2] * diss[ind1] + membr[j2] * diss[ind2] -
                        membr[i2] * membr[j2] * d12 / (membr[i2] + membr[j2])) /
                       (membr[i2] + membr[j2]);
        }
        if (i2 < k) {
          if (diss[ind1] < dmin) {
            dmin = diss[ind1];
            jj = k;
          }
        } else { // i2 > k
          // FIX: the rest of the else clause is a fix by JB to ensure correct
          // nearest neighbours are foudn when a non-monotone clustering
          // method (e.g. the centroid methods) are used
          if (diss[ind1] < disnn[k]) { // find nearest neighbour of i2
            disnn[k] = diss[ind1];
            nn[k] = i2;
          }
        }
      }
    }
    membr[i2] = membr[i2] + membr[j2];
    disnn[i2] = dmin;
    nn[i2] = jj;

    // Update list of NNs insofar as this is required
    for (i = 0; i < n - 1; i++) {
      if (flag[i] && (nn[i] == i2 || nn[i] == j2)) {
        // redetermine nn of i:
        dmin = DBL_MAX;
        for (j = i + 1; j < n; j++) {
          if (flag[j]) {
            ind = ioffst(n, i, j);
            if (diss[ind] < dmin) {
              dmin = diss[ind];
              jj = j;
            }
          }
        }
        // FIXME: if dmin == DBL_MAX then run as fast as you can!
        nn[i] = jj;
        disnn[i] = dmin;
      }
    }
  }
  return 0;
}


/*

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, prepare the seq. of aggloms. and "horiz."    C
C  order of objects for plotting the dendrogram using S routine C
C  'plclust'.                                                   C
C                                                               C
C  Parameters:                                                  C
C                                                               C
C  IA, IB:       vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  IIA, IIB:     used to store IA and IB values differently     C
C                (in form needed for S command 'plclust'        C
C  IORDER:       "horiz." order of objects for dendrogram       C
C                                                               C
C  F. Murtagh, ESA/ESO/STECF, Garching, June 1991               C
C                                                               C
C  HISTORY                                                      C
C                                                               C
C  Adapted from routine HCASS, which additionally determines    C
C   cluster assignments at all levels, at extra comput. expense C
C                                                               C
C  Adapted from fortran to C by Sergio Oller 2022.              C
C  ia and ib inputs are assumed to be 0-based indexed and       C
C  iia and iib are 1-based indexed and use signs as plclust     C
C  expects                                                      C
C---------------------------------------------------------------C
*/

int hcass2(idx_t n, const idx_t *ia, const idx_t *ib, idx_t *iorder, idx_t *iia, idx_t *iib) {
   idx_t i, j, k, k1, k2, loc;
   /*
C     Following bit is to get seq. of merges into format acceptable to plclust
C     I coded clusters as lowest seq. no. of constituents; S's 'hclust' codes
C     singletons as -ve numbers, and non-singletons with their seq. nos.
   */
   for (i=0;i<n;++i) {
      iia[i] = ia[i]+1;
      iib[i] = ib[i]+1;
   }
   for (i=0;i<n-1;++i) {
      //In the following, smallest (+ve or -ve) seq. no. wanted
      k = min(ia[i], ib[i]);
      for (j=i+1;j<n-1;++j) {
         if (ia[j] == k) iia[j] = -(i+1);
         if (ib[j] == k) iib[j] = -(i+1);
      }
   }
   for (i=0;i<n-1;i++) {
      iia[i] *= -1;
      iib[i] *= -1;
   }
   /* Ensure pairs follow:
     - If we have a pair of non-singleton,singleton, permute it so it becomes singleton,non-singleton
     - If both are non singletons use (min,max)
   */
   for (i=0;i<n-1;i++) {
      if (iia[i] > 0 && iib[i] < 0) {
         k = iia[i];
         iia[i] = iib[i];
         iib[i] = k;
      }
      if (iia[i] > 0 && iib[i] > 0) {
         k1 = min(iia[i], iib[i]);
         k2 = max(iia[i], iib[i]);
         iia[i] = k1;
         iib[i] = k2;
      }
   }
   // 'ORDER'
   iorder[0] = iia[n-2];
   iorder[1] = iib[n-2];

   loc=1;
   int found;
   for (i=n-3; i>=0;i--) {
      found = 0;
      for (j=0;j<loc;++j) {
         if (iorder[j] == i+1) {
            // replace iorder[j] with iia[i] and iib[i]
            iorder[j] = iia[i];
            if (j == loc) {
               loc++;
               iorder[loc] = iib[i];
            } else {
               loc++;
               for (k=loc-1;k>j;k--) {
                  iorder[k] = iorder[k-1];
               }
               iorder[j+1] = iib[i];
            }
            found = 1;
            break;
         }
      }
      assert(found != 0);
   }
   for (i=0;i<n;i++) {
      iorder[i] = -iorder[i];
   }
   return 0;
}
