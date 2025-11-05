/*
  This file is part of JURASSIC.
  
  JURASSIC is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  JURASSIC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with JURASSIC. If not, see <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2003-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  JURASSIC retrieval processor.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Carry out optimal estimation retrieval. */
void optimal_estimation(
  ret_t * ret,
  ctl_t * ctl,
  tbl_t * tbl,
  obs_t * obs_meas,
  obs_t * obs_i,
  atm_t * atm_apr,
  atm_t * atm_i);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static atm_t atm_i, atm_apr;
  static ctl_t ctl;
  static obs_t obs_i, obs_meas;
  static ret_t ret;

  FILE *dirlist;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <dirlist>");

  /* Measure CPU-time... */
  TIMER("total", 1);

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  read_ret(argc, argv, &ctl, &ret);

  /* Initialize look-up tables... */
  tbl_t *tbl = read_tbl(&ctl);

  /* Open directory list... */
  if (!(dirlist = fopen(argv[2], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%4999s", ret.dir) != EOF) {

    /* Write info... */
    LOG(1, "\nRetrieve in directory %s...\n", ret.dir);

    /* Read atmospheric data... */
    read_atm(ret.dir, "atm_apr.tab", &ctl, &atm_apr);

    /* Read observation data... */
    read_obs(ret.dir, "obs_meas.tab", &ctl, &obs_meas);

    /* Run retrieval... */
    optimal_estimation(&ret, &ctl, tbl, &obs_meas, &obs_i, &atm_apr, &atm_i);

    /* Measure CPU-time... */
    TIMER("total", 2);
  }

  /* Write info... */
  LOG(1, "\nRetrieval done...");

  /* Measure CPU-time... */
  TIMER("total", 3);

  /* Free... */
  free(tbl);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void optimal_estimation(
  ret_t *ret,
  ctl_t *ctl,
  tbl_t *tbl,
  obs_t *obs_meas,
  obs_t *obs_i,
  atm_t *atm_apr,
  atm_t *atm_i) {

  static int ipa[N], iqa[N];

  FILE *out;

  char filename[2 * LEN];

  double chisq, disq = 0, lmpar = 0.001;

  int it = 0;

  /* ------------------------------------------------------------
     Initialize...
     ------------------------------------------------------------ */

  /* Get sizes... */
  const size_t m = obs2y(ctl, obs_meas, NULL, NULL, NULL);
  const size_t n = atm2x(ctl, atm_apr, NULL, iqa, ipa);
  if (m == 0 || n == 0)
    ERRMSG("Check problem definition!");

  /* Write info... */
  LOG(1, "Problem size: m= %d / n= %d "
      "(alloc= %.4g MB / stat= %.4g MB)",
      (int) m, (int) n,
      (double) (3 * m * n + 4 * n * n + 8 * m +
		8 * n) * sizeof(double) / 1024. / 1024.,
      (double) (5 * sizeof(atm_t) + 3 * sizeof(obs_t)
		+ 2 * N * sizeof(int)) / 1024. / 1024.);

  /* Allocate... */
  gsl_matrix *a = gsl_matrix_alloc(n, n);
  gsl_matrix *cov = gsl_matrix_alloc(n, n);
  gsl_matrix *k_i = gsl_matrix_alloc(m, n);
  gsl_matrix *s_a_inv = gsl_matrix_alloc(n, n);

  gsl_vector *b = gsl_vector_alloc(n);
  gsl_vector *dx = gsl_vector_alloc(n);
  gsl_vector *dy = gsl_vector_alloc(m);
  gsl_vector *sig_eps_inv = gsl_vector_alloc(m);
  gsl_vector *sig_formod = gsl_vector_alloc(m);
  gsl_vector *sig_noise = gsl_vector_alloc(m);
  gsl_vector *x_a = gsl_vector_alloc(n);
  gsl_vector *x_i = gsl_vector_alloc(n);
  gsl_vector *x_step = gsl_vector_alloc(n);
  gsl_vector *y_aux = gsl_vector_alloc(m);
  gsl_vector *y_i = gsl_vector_alloc(m);
  gsl_vector *y_m = gsl_vector_alloc(m);

  /* Set initial state... */
  copy_atm(ctl, atm_i, atm_apr, 0);
  copy_obs(ctl, obs_i, obs_meas, 0);
  formod(ctl, tbl, atm_i, obs_i);

  /* Set state vectors and observation vectors... */
  atm2x(ctl, atm_apr, x_a, NULL, NULL);
  atm2x(ctl, atm_i, x_i, NULL, NULL);
  obs2y(ctl, obs_meas, y_m, NULL, NULL);
  obs2y(ctl, obs_i, y_i, NULL, NULL);

  /* Set inverse a priori covariance S_a^-1... */
  set_cov_apr(ret, ctl, atm_apr, iqa, ipa, s_a_inv);
  write_matrix(ret->dir, "matrix_cov_apr.tab", ctl, s_a_inv,
	       atm_i, obs_i, "x", "x", "r");
  matrix_invert(s_a_inv);

  /* Get measurement errors... */
  set_cov_meas(ret, ctl, obs_meas, sig_noise, sig_formod, sig_eps_inv);

  /* Create cost function file... */
  sprintf(filename, "%s/costs.tab", ret->dir);
  if (!(out = fopen(filename, "w")))
    ERRMSG("Cannot create cost function file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = iteration number\n"
	  "# $2 = normalized cost function\n"
	  "# $3 = number of measurements\n"
	  "# $4 = number of state vector elements\n\n");

  /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
  gsl_vector_memcpy(dx, x_i);
  gsl_vector_sub(dx, x_a);
  gsl_vector_memcpy(dy, y_m);
  gsl_vector_sub(dy, y_i);

  /* Compute cost function... */
  chisq = cost_function(dx, dy, s_a_inv, sig_eps_inv);

  /* Write info... */
  LOG(1, "it= %d / chi^2/m= %g", it, chisq);

  /* Write to cost function file... */
  fprintf(out, "%d %g %d %d\n", it, chisq, (int) m, (int) n);

  /* Compute initial kernel... */
  kernel(ctl, tbl, atm_i, obs_i, k_i);

  /* ------------------------------------------------------------
     Levenberg-Marquardt minimization...
     ------------------------------------------------------------ */

  /* Outer loop... */
  for (it = 1; it <= ret->conv_itmax; it++) {

    /* Store current cost function value... */
    double chisq_old = chisq;

    /* Compute kernel matrix K_i... */
    if (it > 1 && it % ret->kernel_recomp == 0)
      kernel(ctl, tbl, atm_i, obs_i, k_i);

    /* Compute K_i^T * S_eps^{-1} * K_i ... */
    if (it == 1 || it % ret->kernel_recomp == 0)
      matrix_product(k_i, sig_eps_inv, 1, cov);

    /* Determine b = K_i^T * S_eps^{-1} * dy - S_a^{-1} * dx ... */
    for (size_t i = 0; i < m; i++)
      gsl_vector_set(y_aux, i, gsl_vector_get(dy, i)
		     * POW2(gsl_vector_get(sig_eps_inv, i)));
    gsl_blas_dgemv(CblasTrans, 1.0, k_i, y_aux, 0.0, b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, s_a_inv, dx, 1.0, b);

    /* Inner loop... */
    for (int it2 = 0; it2 < 20; it2++) {

      /* Compute A = (1 + lmpar) * S_a^{-1} + K_i^T * S_eps^{-1} * K_i ... */
      gsl_matrix_memcpy(a, s_a_inv);
      gsl_matrix_scale(a, 1 + lmpar);
      gsl_matrix_add(a, cov);

      /* Solve A * x_step = b by means of Cholesky decomposition... */
      gsl_linalg_cholesky_decomp(a);
      gsl_linalg_cholesky_solve(a, b, x_step);

      /* Update atmospheric state... */
      gsl_vector_add(x_i, x_step);
      copy_atm(ctl, atm_i, atm_apr, 0);
      copy_obs(ctl, obs_i, obs_meas, 0);
      x2atm(ctl, x_i, atm_i);

      /* Check atmospheric state... */
      for (int ip = 0; ip < atm_i->np; ip++) {
	atm_i->p[ip] = MIN(MAX(atm_i->p[ip], 5e-7), 5e4);
	atm_i->t[ip] = MIN(MAX(atm_i->t[ip], 100), 400);
	for (int ig = 0; ig < ctl->ng; ig++)
	  atm_i->q[ig][ip] = MIN(MAX(atm_i->q[ig][ip], 0), 1);
	for (int iw = 0; iw < ctl->nw; iw++)
	  atm_i->k[iw][ip] = MAX(atm_i->k[iw][ip], 0);
      }
      atm_i->clz = MAX(atm_i->clz, 0);
      atm_i->cldz = MAX(atm_i->cldz, 0.1);
      for (int icl = 0; icl < ctl->ncl; icl++)
	atm_i->clk[icl] = MAX(atm_i->clk[icl], 0);
      atm_i->sft = MIN(MAX(atm_i->sft, 100), 400);
      for (int isf = 0; isf < ctl->nsf; isf++)
	atm_i->sfeps[isf] = MIN(MAX(atm_i->sfeps[isf], 0), 1);

      /* Forward calculation... */
      formod(ctl, tbl, atm_i, obs_i);
      obs2y(ctl, obs_i, y_i, NULL, NULL);

      /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
      gsl_vector_memcpy(dx, x_i);
      gsl_vector_sub(dx, x_a);
      gsl_vector_memcpy(dy, y_m);
      gsl_vector_sub(dy, y_i);

      /* Compute cost function... */
      chisq = cost_function(dx, dy, s_a_inv, sig_eps_inv);

      /* Modify Levenberg-Marquardt parameter... */
      if (chisq > chisq_old) {
	lmpar *= 10;
	gsl_vector_sub(x_i, x_step);
      } else {
	lmpar /= 10;
	break;
      }
    }

    /* Write info... */
    LOG(1, "it= %d / chi^2/m= %g", it, chisq);

    /* Write to cost function file... */
    fprintf(out, "%d %g %d %d\n", it, chisq, (int) m, (int) n);

    /* Get normalized step size in state space... */
    gsl_blas_ddot(x_step, b, &disq);
    disq /= (double) n;

    /* Convergence test... */
    if ((it == 1 || it % ret->kernel_recomp == 0) && disq < ret->conv_dmin)
      break;
  }

  /* Close cost function file... */
  fclose(out);

  /* Store results... */
  write_atm(ret->dir, "atm_final.tab", ctl, atm_i);
  write_obs(ret->dir, "obs_final.tab", ctl, obs_i);
  write_matrix(ret->dir, "matrix_kernel.tab", ctl, k_i,
	       atm_i, obs_i, "y", "x", "r");

  /* ------------------------------------------------------------
     Analysis of retrieval results...
     ------------------------------------------------------------ */

  /* Check if error analysis is requested... */
  if (ret->err_ana) {

    /* Allocate... */
    gsl_matrix *auxnm = gsl_matrix_alloc(n, m);
    gsl_matrix *corr = gsl_matrix_alloc(n, n);
    gsl_matrix *gain = gsl_matrix_alloc(n, m);

    /* Compute inverse retrieval covariance...
       cov^{-1} = S_a^{-1} + K_i^T * S_eps^{-1} * K_i */
    matrix_product(k_i, sig_eps_inv, 1, cov);
    gsl_matrix_add(cov, s_a_inv);

    /* Compute retrieval covariance... */
    matrix_invert(cov);
    write_matrix(ret->dir, "matrix_cov_ret.tab", ctl, cov,
		 atm_i, obs_i, "x", "x", "r");
    write_stddev("total", ret, ctl, atm_i, cov);

    /* Compute correlation matrix... */
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	gsl_matrix_set(corr, i, j, gsl_matrix_get(cov, i, j)
		       / sqrt(gsl_matrix_get(cov, i, i))
		       / sqrt(gsl_matrix_get(cov, j, j)));
    write_matrix(ret->dir, "matrix_corr.tab", ctl, corr,
		 atm_i, obs_i, "x", "x", "r");

    /* Compute gain matrix...
       G = cov * K^T * S_eps^{-1} */
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < m; j++)
	gsl_matrix_set(auxnm, i, j, gsl_matrix_get(k_i, j, i)
		       * POW2(gsl_vector_get(sig_eps_inv, j)));
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cov, auxnm, 0.0, gain);
    write_matrix(ret->dir, "matrix_gain.tab", ctl, gain,
		 atm_i, obs_i, "x", "y", "c");

    /* Compute retrieval error due to noise... */
    matrix_product(gain, sig_noise, 2, a);
    write_stddev("noise", ret, ctl, atm_i, a);

    /* Compute retrieval error  due to forward model errors... */
    matrix_product(gain, sig_formod, 2, a);
    write_stddev("formod", ret, ctl, atm_i, a);

    /* Compute averaging kernel matrix
       A = G * K ... */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, gain, k_i, 0.0, a);
    write_matrix(ret->dir, "matrix_avk.tab", ctl, a,
		 atm_i, obs_i, "x", "x", "r");

    /* Analyze averaging kernel matrix... */
    analyze_avk(ret, ctl, atm_i, iqa, ipa, a);

    /* Free... */
    gsl_matrix_free(auxnm);
    gsl_matrix_free(corr);
    gsl_matrix_free(gain);
  }

  /* ------------------------------------------------------------
     Finalize...
     ------------------------------------------------------------ */

  gsl_matrix_free(a);
  gsl_matrix_free(cov);
  gsl_matrix_free(k_i);
  gsl_matrix_free(s_a_inv);

  gsl_vector_free(b);
  gsl_vector_free(dx);
  gsl_vector_free(dy);
  gsl_vector_free(sig_eps_inv);
  gsl_vector_free(sig_formod);
  gsl_vector_free(sig_noise);
  gsl_vector_free(x_a);
  gsl_vector_free(x_i);
  gsl_vector_free(x_step);
  gsl_vector_free(y_aux);
  gsl_vector_free(y_i);
  gsl_vector_free(y_m);
}
