/*
 * check_baseline.h
 *
 * The bl_trespass is used to save the details about the
 * baseline and patch exceptions. When a particular baseline
 * happens to stick outside the borders of the visibility
 * uv grid, or a particular patch is too large to fit into 
 * the UV grid, the information about the station pair and 
 * the frequency is stored in bl_trespass record.
 * Tre records make up a table with N_TRESPASSED_BASELINES
 * lines, which is later used to generate a diagnostic report.
 *
 * Created by L. Benkevitch, 03-Dec-2010 
 */

#define N_TRESPASSED_BASELINES 1024

struct bl_trespass
{
  int   st1;      /* Station 1 */
  int   st2;      /* Station 2 */
  double u_len;   /* U in meters */
  double v_len;   /* V in meters */
  double freq;    /* (Minimum) frequency */
  double gha;     /* Greenwich Hour Angle (GHA) in radians */
  int   ch;       /* Frequency channel (IF) */
  int   cch;      /* Correlator freq channels within an IF. */
  int   usize;    /* Patch size in cells along U  */
  int   vsize;    /* Patch size in cells along V */ 
  int   ercode;   /* Error code as given in get_patch() */
};
