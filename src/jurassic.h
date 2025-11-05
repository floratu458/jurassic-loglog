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
  JURASSIC library declarations.
*/

/*! 
  \mainpage
  
  The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast
  infrared radiative transfer model for the analysis of atmospheric
  remote sensing measurements.

  \section Introduction

  The source code of JURASSIC is available from the
  [git repository](https://github.com/slcs-jsc/jurassic). Please see the
  [README.md](https://github.com/slcs-jsc/jurassic/blob/master/README.md)
  in the git repository for introductory information. More information
  can be found in the [user manual](https://slcs-jsc.github.io/jurassic).
  
  This doxygen manual contains information about the algorithms and
  data structures used in the code. Please refer to the `jurassic.h'
  documentation for a first overview.
  
  \section References
  
  For citing the model in scientific publications, please see
  [CITATION.cff](https://github.com/slcs-jsc/jurassic/blob/master/CITATION.cff)
  and refer to the following papers:
  
  _Baumeister, P. F. and Hoffmann, L.: Fast infrared radiative
  transfer calculations using graphics processing units: JURASSIC-GPU
  v2.0, Geosci. Model Dev., 15, 1855–1874,
  https://doi.org/10.5194/gmd-15-1855-2022, 2022._
  
  _Hoffmann, L., and M. J. Alexander, Retrieval of stratospheric
  temperatures from Atmospheric Infrared Sounder radiance measurements
  for gravity wave studies, J. Geophys. Res., 114, D07105,
  https://doi.org/10.1029/2008JD011241, 2009._
  
  _Hoffmann, L., Kaufmann, M., Spang, R., Müller, R., Remedios, J. J.,
  Moore, D. P., Volk, C. M., von Clarmann, T., and Riese, M.: Envisat
  MIPAS measurements of CFC-11: retrieval, validation, and
  climatology, Atmos. Chem. Phys., 8, 3671-3688,
  https://doi.org/10.5194/acp-8-3671-2008, 2008._
  
  Additional references are collected here:
  https://slcs-jsc.github.io/jurassic/references
  
  \section License
  
  JURASSIC is being develop at the Jülich Supercomputing Centre,
  Forschungszentrum Jülich, Germany.
  
  JURASSIC is distributed under the terms of the
  [GNU General Public License v3.0](https://github.com/slcs-jsc/jurassic/blob/master/COPYING).
  
  \section Contributing
  
  We are interested in supporting operational and research
  applications with JURASSIC.
  
  You can submit bug reports or feature requests on the
  [issue tracker](https://github.com/slcs-jsc/jurassic/issues).
  
  Proposed code changes and fixes can be submitted as
  [pull requests](https://github.com/slcs-jsc/jurassic/pulls).
  
  Please do not hesitate to contact us if you have any questions or
  need assistance.
  
  \section Contact
  
  Dr. Lars Hoffmann
  
  Jülich Supercomputing Centre, Forschungszentrum Jülich
  
  e-mail: <l.hoffmann@fz-juelich.de>
*/

#ifndef JURASSIC_H
#define JURASSIC_H

/* ------------------------------------------------------------
   Includes...
   ------------------------------------------------------------ */

#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#ifndef C1
#define C1 1.19104259e-8
#endif

/*! Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#ifndef C2
#define C2 1.43877506
#endif

/*! Minimum emissivity. */
#ifndef EPSMIN
#define EPSMIN 0
#endif

/*! Maximum emissivity. */
#ifndef EPSMAX
#define EPSMAX 1
#endif

/*! Standard gravity [m/s^2]. */
#ifndef G0
#define G0 9.80665
#endif

/*! Standard scale height [km]. */
#ifndef H0
#define H0 7.0
#endif

/*! Boltzmann constant [kg m^2/(K s^2)]. */
#ifndef KB
#define KB 1.3806504e-23
#endif

/*! Mass of Earth [kg]. */
#ifndef ME
#define ME 5.976e24
#endif

/*! Avogadro's number. */
#ifndef NA
#define NA 6.02214199e23
#endif

/*! Nitrogen concentration. */
#ifndef N2
#define N2 0.78084
#endif

/*! Oxygen concentration. */
#ifndef O2
#define O2 0.20946
#endif

/*! Standard pressure [hPa]. */
#ifndef P0
#define P0 1013.25
#endif

/*! Mean radius of Earth [km]. */
#ifndef RE
#define RE 6367.421
#endif

/*! Ideal gas constant [J/(mol K)]. */
#ifndef RI
#define RI 8.3144598
#endif

/*! Standard temperature [K]. */
#ifndef T0
#define T0 273.15
#endif

/*! Minimum temperature for source function [K]. */
#ifndef TMIN
#define TMIN 100.
#endif

/*! Maximum temperature for source function [K]. */
#ifndef TMAX
#define TMAX 400.
#endif

/*! Effective temperature of the sun [K]. */
#ifndef TSUN
#define TSUN 5780.
#endif

/*! Minimum column density [molecules/cm^2]. */
#ifndef UMIN
#define UMIN 0
#endif

/*! Maximum column density [molecules/cm^2]. */
#ifndef UMAX
#define UMAX 1e30
#endif

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of cloud layer spectral grid points. */
#ifndef NCL
#define NCL 8
#endif

/*! Maximum number of radiance channels. */
#ifndef ND
#define ND 128
#endif

/*! Maximum number of emitters. */
#ifndef NG
#define NG 8
#endif

/*! Maximum number of atmospheric data points. */
#ifndef NP
#define NP 256
#endif

/*! Maximum number of ray paths. */
#ifndef NR
#define NR 256
#endif

/*! Maximum number of surface layer spectral grid points. */
#ifndef NSF
#define NSF 8
#endif

/*! Maximum number of spectral windows. */
#ifndef NW
#define NW 4
#endif

/*! Maximum length of ASCII data lines. */
#ifndef LEN
#define LEN 10000
#endif

/*! Maximum size of measurement vector. */
#ifndef M
#define M (NR*ND)
#endif

/*! Maximum size of state vector. */
#ifndef N
#define N ((2 + NG + NW) * NP + NCL + NSF + 3)
#endif

/*! Maximum number of quantities. */
#ifndef NQ
#define NQ (5 + NG + NW + NCL + NSF)
#endif

/*! Maximum number of LOS points. */
#ifndef NLOS
#define NLOS 4096
#endif

/*! Maximum number of shape function grid points. */
#ifndef NSHAPE
#define NSHAPE 20000
#endif

/*! Number of ray paths used for FOV calculations. */
#ifndef NFOV
#define NFOV 5
#endif

/*! Maximum number of pressure levels in emissivity tables. */
#ifndef TBLNP
#define TBLNP 41
#endif

/*! Maximum number of temperatures in emissivity tables. */
#ifndef TBLNT
#define TBLNT 30
#endif

/*! Maximum number of column densities in emissivity tables. */
#ifndef TBLNU
#define TBLNU 320
#endif

/*! Maximum number of source function temperature levels. */
#ifndef TBLNS
#define TBLNS 1200
#endif

/*! Maximum number of RFM spectral grid points. */
#ifndef RFMNPTS
#define RFMNPTS 10000000
#endif

/* ------------------------------------------------------------
   Quantity indices...
   ------------------------------------------------------------ */

/*! Index for pressure. */
#define IDXP 0

/*! Index for temperature. */
#define IDXT 1

/*! Indices for volume mixing ratios. */
#define IDXQ(ig) (2 + (ig))

/*! Indices for extinction. */
#define IDXK(iw) (2 + (ctl->ng) + (iw))

/*! Index for cloud layer height. */
#define IDXCLZ (2 + (ctl->ng) + (ctl->nw))

/*! Index for cloud layer depth. */
#define IDXCLDZ (3 + (ctl->ng) + (ctl->nw))

/*! Indices for cloud layer extinction. */
#define IDXCLK(icl) (4 + (ctl->ng) + (ctl->nw) + (icl))

/*! Index for surface layer temperature. */
#define IDXSFT (4 + (ctl->ng) + (ctl->nw) + (ctl->ncl))

/*! Indices for surface layer emissivity. */
#define IDXSFEPS(isf) (5 + (ctl->ng) + (ctl->nw) + (ctl->ncl) + (isf))

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/**
 * @brief Allocate memory for an array.
 *
 * Allocates a contiguous block of memory for an array of `n` elements of
 * type `type` and assigns the pointer to `ptr`. If allocation fails, the
 * program prints an error message and terminates.
 *
 * @param[out] ptr Pointer to be allocated.
 * @param[in] type Data type of the elements to allocate.
 * @param[in] n Number of elements to allocate.
 *
 * @note Wraps `malloc()` with built-in error handling.
 *
 * @see FREAD, FWRITE
 *
 * @author Lars Hoffmann
 */
#define ALLOC(ptr, type, n) \
  if((ptr=malloc((size_t)(n)*sizeof(type)))==NULL) \
    ERRMSG("Out of memory!");

/**
 * @brief Compute brightness temperature from radiance.
 *
 * Computes the equivalent blackbody (brightness) temperature corresponding to
 * a given spectral radiance and wavenumber, using the inverse Planck function.
 * This form assumes the spectroscopic constants `C1` and `C2` are defined for
 * wavenumber units (cm⁻¹).
 *
 * @param[in] rad Spectral radiance [W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹].
 * @param[in] nu  Wavenumber [cm⁻¹].
 *
 * @return Brightness temperature [K].
 *
 * @see PLANCK, C1, C2
 *
 * @note Based on Planck’s law in wavenumber form:
 *       \f$ T_b = \frac{c_2 \nu}{\ln(1 + \frac{c_1 \nu^3}{L_\nu})} \f$
 *       where \( L_\nu \) is radiance.
 *
 * @author Lars Hoffmann
 */
#define BRIGHT(rad, nu) \
  (C2 * (nu) / gsl_log1p(C1 * POW3(nu) / (rad)))

/**
 * @brief Convert degrees to radians.
 *
 * Converts an angle measured in degrees to radians using:
 * \f$ \text{radians} = \text{degrees} \times \frac{\pi}{180} \f$.
 *
 * @param[in] deg Angle in degrees.
 *
 * @return Angle in radians.
 *
 * @see RAD2DEG
 *
 * @author Lars Hoffmann
 */
#define DEG2RAD(deg) ((deg) * (M_PI / 180.0))

/**
 * @brief Compute Cartesian distance between two 3D vectors.
 *
 * Computes the Euclidean distance between two 3D vectors or points.
 *
 * @param[in] a First vector (array of length 3).
 * @param[in] b Second vector (array of length 3).
 *
 * @return Euclidean distance between a and b.
 *
 * @see DIST2
 *
 * @note Equivalent to \f$ \sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2} \f$.
 *
 * @author Lars Hoffmann
 */
#define DIST(a, b) sqrt(DIST2(a, b))

/**
 * @brief Compute squared distance between two 3D vectors.
 *
 * Computes the square of the Euclidean distance between two 3D vectors.
 * Useful when only relative distances are needed (avoids square root).
 *
 * @param[in] a First vector (array of length 3).
 * @param[in] b Second vector (array of length 3).
 *
 * @return Squared distance between a and b.
 *
 * @see DIST
 *
 * @author Lars Hoffmann
 */
#define DIST2(a, b) \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/**
 * @brief Compute dot product of two 3D vectors.
 *
 * Computes the scalar (dot) product between two 3-element vectors:
 * \f$ a \cdot b = a_x b_x + a_y b_y + a_z b_z \f$.
 *
 * @param[in] a First vector (array of length 3).
 * @param[in] b Second vector (array of length 3).
 *
 * @return Scalar dot product of a and b.
 *
 * @see NORM
 *
 * @author Lars Hoffmann
 */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/**
 * @brief Read binary data from a file.
 *
 * Reads `size` elements of type `type` from the file stream `out` into
 * the memory pointed to by `ptr`. If the number of elements read differs
 * from `size`, an error message is printed and the program terminates.
 *
 * @param[out] ptr Pointer to destination memory buffer.
 * @param[in] type Type of each element to read.
 * @param[in] size Number of elements to read.
 * @param[in,out] out File stream opened for reading.
 *
 * @see FWRITE, ALLOC
 *
 * @note Wraps `fread()` with error handling.
 *
 * @author Lars Hoffmann
 */
#define FREAD(ptr, type, size, out) { \
    if(fread(ptr, sizeof(type), size, out)!=size) \
      ERRMSG("Error while reading!"); \
  }

/**
 * @brief Write binary data to a file.
 *
 * Writes `size` elements of type `type` from the memory pointed to by `ptr`
 * to the file stream `out`. If the number of elements written differs
 * from `size`, an error message is printed and the program terminates.
 *
 * @param[in] ptr Pointer to memory buffer containing data to write.
 * @param[in] type Type of each element to write.
 * @param[in] size Number of elements to write.
 * @param[in,out] out File stream opened for writing.
 *
 * @see FREAD
 *
 * @note Wraps `fwrite()` with error handling.
 *
 * @author Lars Hoffmann
 */
#define FWRITE(ptr, type, size, out) { \
    if(fwrite(ptr, sizeof(type), size, out)!=size) \
      ERRMSG("Error while writing!"); \
  }

/**
 * @brief Determine the maximum of two values.
 *
 * Returns the greater of two scalar values.
 *
 * @param[in] a First value.
 * @param[in] b Second value.
 *
 * @return Maximum of a and b.
 *
 * @see MIN
 *
 * @note Both arguments are evaluated multiple times; avoid side effects.
 *
 * @author Lars Hoffmann
 */
#define MAX(a,b) (((a)>(b))?(a):(b))

/**
 * @brief Determine the minimum of two values.
 *
 * Returns the smaller of two scalar values.
 *
 * @param[in] a First value.
 * @param[in] b Second value.
 *
 * @return Minimum of a and b.
 *
 * @see MAX
 *
 * @note Both arguments are evaluated multiple times; avoid side effects.
 *
 * @author Lars Hoffmann
 */
#define MIN(a,b) (((a)<(b))?(a):(b))

/**
 * @brief Compute linear interpolation.
 *
 * Performs simple linear interpolation between two known points
 * (x₀, y₀) and (x₁, y₁) to estimate the value of y at a given x.
 *
 * @param[in] x0 Lower x-value.
 * @param[in] y0 Function value at x₀.
 * @param[in] x1 Upper x-value.
 * @param[in] y1 Function value at x₁.
 * @param[in] x Interpolation point.
 *
 * @return Linearly interpolated y-value at x.
 *
 * @see LOGX, LOGY
 *
 * @author Lars Hoffmann
 */
#define LIN(x0, y0, x1, y1, x) \
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/**
 * @brief Compute logarithmic interpolation in x.
 *
 * Performs interpolation assuming logarithmic variation in the x-axis.
 * If either x/x₀ or x₁/x₀ is nonpositive, reverts to linear interpolation.
 *
 * @param[in] x0 Lower x-value.
 * @param[in] y0 Function value at x₀.
 * @param[in] x1 Upper x-value.
 * @param[in] y1 Function value at x₁.
 * @param[in] x Interpolation point.
 *
 * @return Interpolated y-value at x.
 *
 * @see LIN, LOGY
 *
 * @author Lars Hoffmann
 */
#define LOGX(x0, y0, x1, y1, x) \
  (((x)/(x0)>0 && (x1)/(x0)>0) \
   ? ((y0)+((y1)-(y0))*log((x)/(x0))/log((x1)/(x0))) \
   : LIN(x0, y0, x1, y1, x))

/**
 * @brief Compute logarithmic interpolation in y.
 *
 * Performs interpolation assuming exponential variation in the y-axis
 * (logarithmic in y). If y₁/y₀ is nonpositive, reverts to linear interpolation.
 *
 * @param[in] x0 Lower x-value.
 * @param[in] y0 Function value at x₀.
 * @param[in] x1 Upper x-value.
 * @param[in] y1 Function value at x₁.
 * @param[in] x Interpolation point.
 *
 * @return Interpolated y-value at x.
 *
 * @see LIN, LOGX
 *
 * @author Lars Hoffmann
 */
#define LOGY(x0, y0, x1, y1, x) \
  (((y1)/(y0)>0) \
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0)))) \
   : LIN(x0, y0, x1, y1, x))

/**
 * @brief Compute the norm (magnitude) of a 3D vector.
 *
 * Computes the Euclidean norm using the dot product:
 * \f$ |a| = \sqrt{a \cdot a} \f$.
 *
 * @param[in] a Input vector (array of length 3).
 *
 * @return Magnitude (norm) of vector a.
 *
 * @see DOTP
 *
 * @author Lars Hoffmann
 */
#define NORM(a) sqrt(DOTP(a, a))

/**
 * @brief Compute the Planck function in wavenumber form.
 *
 * Computes spectral radiance per unit wavenumber using Planck’s law:
 * \f$ B_\nu(T) = \frac{C_1 \nu^3}{\exp(C_2 \nu / T) - 1} \f$.
 *
 * @param[in] T Temperature [K].
 * @param[in] nu Wavenumber [cm⁻¹].
 *
 * @return Spectral radiance [W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹].
 *
 * @see BRIGHT, C1, C2
 *
 * @note Constants `C1` and `C2` must correspond to wavenumber units.
 *
 * @author Lars Hoffmann
 */
#define PLANCK(T, nu) \
  (C1 * POW3(nu) / gsl_expm1(C2 * (nu) / (T)))

/**
 * @brief Compute the square of a value.
 *
 * Returns x². Inline alternative to `pow(x,2)`.
 *
 * @param[in] x Input value.
 *
 * @return x squared.
 *
 * @author Lars Hoffmann
 */
#define POW2(x) ((x)*(x))

/**
 * @brief Compute the cube of a value.
 *
 * Returns x³. Inline alternative to `pow(x,3)`.
 *
 * @param[in] x Input value.
 *
 * @return x cubed.
 *
 * @author Lars Hoffmann
 */
#define POW3(x) ((x)*(x)*(x))

/**
 * @brief Convert radians to degrees.
 *
 * Converts an angle measured in radians to degrees using:
 * \f$ \text{degrees} = \text{radians} \times \frac{180}{\pi} \f$.
 *
 * @param[in] rad Angle in radians.
 *
 * @return Angle in degrees.
 *
 * @see DEG2RAD
 *
 * @author Lars Hoffmann
 */
#define RAD2DEG(rad) ((rad) * (180.0 / M_PI))

/**
 * @brief Compute air refractivity (n - 1).
 *
 * Approximates the refractivity of air under standard conditions using:
 * \f$ n - 1 = 7.753\times10^{-5} \frac{p}{T} \f$,
 * where p is pressure (hPa) and T is temperature (K).
 *
 * @param[in] p Pressure [hPa].
 * @param[in] T Temperature [K].
 *
 * @return Refractivity (dimensionless, n - 1).
 *
 * @author Lars Hoffmann
 */
#define REFRAC(p, T) (7.753e-05 * (p) / (T))

/**
 * @brief Start or stop a named timer.
 *
 * Calls the `timer()` function with contextual information (file name, function
 * name, and line number) to start or stop timing. Useful for performance profiling.
 *
 * @param[in] name Name or label of the timer.
 * @param[in] mode Operation mode (e.g., start or stop).
 *
 * @note Relies on a user-defined `timer()` function.
 *
 * @author Lars Hoffmann
 */
#define TIMER(name, mode) \
  {timer(name, __FILE__, __func__, __LINE__, mode);}

/**
 * @brief Tokenize a string and parse a variable.
 *
 * Splits a text line into tokens separated by spaces or tabs, and reads
 * a value from the first token using `sscanf()` and the provided format.
 * If tokenization or parsing fails, an error message is printed.
 *
 * @param[in,out] line Input string buffer to tokenize (modified by strtok()).
 * @param[out] tok Pointer to token string.
 * @param[in] format Format string for `sscanf()`.
 * @param[out] var Variable to store parsed value.
 *
 * @note Uses `strtok()` internally and modifies the input buffer.
 *
 * @see FREAD, FWRITE
 *
 * @author Lars Hoffmann
 */
#define TOK(line, tok, format, var) { \
    if(((tok)=strtok((line), " \t"))) { \
      if(sscanf(tok, format, &(var))!=1) continue; \
    } else ERRMSG("Error while reading!"); \
  }

/* ------------------------------------------------------------
   Log messages...
   ------------------------------------------------------------ */

/*! Level of log messages (0=none, 1=basic, 2=detailed, 3=debug). */
#ifndef LOGLEV
#define LOGLEV 2
#endif

/*!
 * \brief Print a log message with a specified logging level.
 *
 * This macro prints a formatted log message to the standard output if
 * the specified logging level meets certain conditions. The message
 * will be indented if the logging level is greater than or equal to
 * 2.
 * 
 * \param level The logging level of the message. This should be an integer value.
 * \param ... The formatted message string and its arguments, similar to printf.
 *
 * \details
 * The `LOG` macro provides a simple way to log messages with
 * different levels of importance. The message is only printed if the
 * specified `level` is less than or equal to the pre-defined `LOGLEV`
 * macro. If the `level` is greater than or equal to 2, the message is
 * preceded by two spaces for indentation.
 *
 * The macro expands to a block of code that:
 * - Checks if the `level` is greater than or equal to 2, and if so, prints two spaces.
 * - Checks if the `level` is less than or equal to `LOGLEV`, and if so, prints the
 *   formatted message followed by a newline.
 *
 * \note
 * The `LOGLEV` macro must be defined with an appropriate logging level
 * before using the `LOG` macro.
 * 
 * @author Lars Hoffmann
 */
#define LOG(level, ...) {						\
    if(level >= 2)							\
      printf("  ");							\
    if(level <= LOGLEV) {						\
      printf(__VA_ARGS__);						\
      printf("\n");							\
    }									\
  }

/*!
 * \brief Print a warning message with contextual information.
 *
 * This macro prints a formatted warning message to the standard
 * output, including the file name, function name, and line number
 * where the warning occurred. The message is then passed to the `LOG`
 * macro with a logging level of 0.
 * 
 * \param ... The formatted warning message string and its arguments, similar to printf.
 *
 * \details
 * The `WARN` macro is used to print warning messages with additional context
 * about where the warning was triggered. The message includes the following
 * contextual information:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 *
 * After printing this contextual information, the macro uses the
 * `LOG` macro with a logging level of 0 to print the actual warning
 * message. This ensures that warning messages are always logged,
 * regardless of the value of `LOGLEV`.
 *
 * \note
 * The `LOG` macro must be defined before using the `WARN` macro.
 * 
 * @author Lars Hoffmann
 */
#define WARN(...) {							\
    printf("\nWarning (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
  }

/*!
 * \brief Print an error message with contextual information and terminate the program.
 *
 * This macro prints a formatted error message to the standard output,
 * including the file name, function name, and line number where the
 * error occurred. After printing the message, the program is
 * terminated with an exit status indicating failure.
 * 
 * \param ... The formatted error message string and its arguments, similar to printf.
 *
 * \details
 * The `ERRMSG` macro is used to report critical errors that require the
 * program to terminate immediately. The message includes the following
 * contextual information:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 *
 * After printing this contextual information, the macro uses the
 * `LOG` macro with a logging level of 0 to print the actual error
 * message. Finally, the program exits with a failure status
 * (`EXIT_FAILURE`).
 *
 * \note
 * The `LOG` macro must be defined before using the `ERRMSG` macro.
 * 
 * @author Lars Hoffmann
 */
#define ERRMSG(...) {							\
    printf("\nError (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
    exit(EXIT_FAILURE);							\
  }

/*!
 * \brief Print the value of a variable with contextual information.
 *
 * This macro prints the value of a variable to the standard output,
 * including the file name, function name, and line number where the
 * macro is called. The output also includes the variable's name and
 * value in a formatted string.
 * 
 * \param format The format string used to print the variable's value, similar to printf.
 * \param var The variable to be printed.
 *
 * \details
 * The `PRINT` macro is used to output the value of a variable along with
 * additional context about where the macro is called. The message includes:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 * - The name of the variable being printed (`#var`).
 * - The value of the variable, formatted according to the provided format string (`format`).
 *
 * This macro is particularly useful for debugging purposes, providing
 * a convenient way to trace variable values and their locations in
 * the code.
 *
 * \note
 * The format string must be compatible with the type of the variable being printed.
 * 
 * @author Lars Hoffmann
 */
#define PRINT(format, var)						\
  printf("Print (%s, %s, l%d): %s= "format"\n",				\
	 __FILE__, __func__, __LINE__, #var, var);

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/**
 * @brief Atmospheric profile data.
 *
 * Holds one vertical atmospheric column including geolocation,
 * thermodynamic, cloud, and surface properties for radiative-transfer
 * calculations.
 */
typedef struct {

  /*! Number of data points. */
  int np;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NP];

  /*! Altitude [km]. */
  double z[NP];

  /*! Longitude [deg]. */
  double lon[NP];

  /*! Latitude [deg]. */
  double lat[NP];

  /*! Pressure [hPa]. */
  double p[NP];

  /*! Temperature [K]. */
  double t[NP];

  /*! Volume mixing ratio [ppv]. */
  double q[NG][NP];

  /*! Extinction [km^-1]. */
  double k[NW][NP];

  /*! Cloud layer height [km]. */
  double clz;

  /*! Cloud layer depth [km]. */
  double cldz;

  /*! Cloud layer extinction [km^-1]. */
  double clk[NCL];

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[NSF];

} atm_t;

/**
 * @brief Control parameters.
 * 
 * This structure contains all control parameters used by the JURASSIC
 * model. The struct is used to collect and to easily pass the control
 * parameters on to the various functions.
 */
typedef struct {

  /*! Number of emitters. */
  int ng;

  /*! Name of each emitter. */
  char emitter[NG][LEN];

  /*! Emitter index of CO2. */
  int ig_co2;

  /*! Emitter index of H2O. */
  int ig_h2o;

  /*! Emitter index of N2. */
  int ig_n2;

  /*! Emitter index of O2. */
  int ig_o2;

  /*! Number of radiance channels. */
  int nd;

  /*! Centroid wavenumber of each channel [cm^-1]. */
  double nu[ND];

  /*! Number of spectral windows. */
  int nw;

  /*! Window index of each channel. */
  int window[ND];

  /*! Number of cloud layer spectral grid points. */
  int ncl;

  /*! Cloud layer wavenumber [cm^-1]. */
  double clnu[NCL];

  /*! Number of surface layer spectral grid points. */
  int nsf;

  /*! Surface layer wavenumber [cm^-1]. */
  double sfnu[NSF];

  /*! Surface treatment (0=none, 1=emissions, 2=downward, 3=solar). */
  int sftype;

  /*! Solar zenith angle at the surface [deg] (-999=auto). */
  double sfsza;

  /*! Basename for table files and filter function files. */
  char tblbase[LEN];

  /*! Look-up table file format (1=ASCII, 2=binary). */
  int tblfmt;

  /*! Reference height for hydrostatic pressure profile (-999 to skip) [km]. */
  double hydz;

  /*! Compute CO2 continuum (0=no, 1=yes). */
  int ctm_co2;

  /*! Compute H2O continuum (0=no, 1=yes). */
  int ctm_h2o;

  /*! Compute N2 continuum (0=no, 1=yes). */
  int ctm_n2;

  /*! Compute O2 continuum (0=no, 1=yes). */
  int ctm_o2;

  /*! Take into account refractivity (0=no, 1=yes). */
  int refrac;

  /*! Maximum step length for raytracing [km]. */
  double rayds;

  /*! Vertical step length for raytracing [km]. */
  double raydz;

  /*! Field-of-view data file. */
  char fov[LEN];

  /*! Field-of-view vertical distance [km]. */
  double fov_dz[NSHAPE];

  /*! Field-of-view weighting factor. */
  double fov_w[NSHAPE];

  /*! Field-of-view number of data points. */
  int fov_n;

  /*! Minimum altitude for pressure retrieval [km]. */
  double retp_zmin;

  /*! Maximum altitude for pressure retrieval [km]. */
  double retp_zmax;

  /*! Minimum altitude for temperature retrieval [km]. */
  double rett_zmin;

  /*! Maximum altitude for temperature retrieval [km]. */
  double rett_zmax;

  /*! Minimum altitude for volume mixing ratio retrieval [km]. */
  double retq_zmin[NG];

  /*! Maximum altitude for volume mixing ratio retrieval [km]. */
  double retq_zmax[NG];

  /*! Minimum altitude for extinction retrieval [km]. */
  double retk_zmin[NW];

  /*! Maximum altitude for extinction retrieval [km]. */
  double retk_zmax[NW];

  /*! Retrieve cloud layer height (0=no, 1=yes). */
  int ret_clz;

  /*! Retrieve cloud layer depth (0=no, 1=yes). */
  int ret_cldz;

  /*! Retrieve cloud layer extinction (0=no, 1=yes). */
  int ret_clk;

  /*! Retrieve surface layer temperature (0=no, 1=yes). */
  int ret_sft;

  /*! Retrieve surface layer emissivity (0=no, 1=yes). */
  int ret_sfeps;

  /*! Use brightness temperature instead of radiance (0=no, 1=yes). */
  int write_bbt;

  /*! Write matrix file (0=no, 1=yes). */
  int write_matrix;

  /*! Forward model (0=CGA, 1=EGA, 2=RFM). */
  int formod;

  /*! Path to RFM binary. */
  char rfmbin[LEN];

  /*! HITRAN file for RFM. */
  char rfmhit[LEN];

  /*! Emitter cross-section files for RFM. */
  char rfmxsc[NG][LEN];

} ctl_t;

/**
 * @brief Line-of-sight data.
 *
 * Contains all quantities along a ray path used for radiative-transfer
 * calculations, including geometry, thermodynamic state, gas and
 * extinction profiles, and precomputed optical parameters.
 */
typedef struct {

  /*! Number of LOS points. */
  int np;

  /*! Altitude [km]. */
  double z[NLOS];

  /*! Longitude [deg]. */
  double lon[NLOS];

  /*! Latitude [deg]. */
  double lat[NLOS];

  /*! Pressure [hPa]. */
  double p[NLOS];

  /*! Temperature [K]. */
  double t[NLOS];

  /*! Volume mixing ratio [ppv]. */
  double q[NLOS][NG];

  /*! Extinction [km^-1]. */
  double k[NLOS][ND];

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[ND];

  /*! Segment length [km]. */
  double ds[NLOS];

  /*! Column density [molecules/cm^2]. */
  double u[NLOS][NG];

  /*! Curtis-Godson pressure [hPa]. */
  double cgp[NLOS][NG];

  /*! Curtis-Godson temperature [K]. */
  double cgt[NLOS][NG];

  /*! Curtis-Godson column density [molecules/cm^2]. */
  double cgu[NLOS][NG];

  /*! Segment emissivity. */
  double eps[NLOS][ND];

  /*! Segment source function [W/(m^2 sr cm^-1)]. */
  double src[NLOS][ND];

} los_t;

/**
 * @brief Observation geometry and radiance data.
 *
 * Stores viewing geometry and radiative quantities for multiple ray paths.
 * Each path represents a line of sight between observer and tangent point,
 * including associated time and location data.
 */
typedef struct {

  /*! Number of ray paths. */
  int nr;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NR];

  /*! Observer altitude [km]. */
  double obsz[NR];

  /*! Observer longitude [deg]. */
  double obslon[NR];

  /*! Observer latitude [deg]. */
  double obslat[NR];

  /*! View point altitude [km]. */
  double vpz[NR];

  /*! View point longitude [deg]. */
  double vplon[NR];

  /*! View point latitude [deg]. */
  double vplat[NR];

  /*! Tangent point altitude [km]. */
  double tpz[NR];

  /*! Tangent point longitude [deg]. */
  double tplon[NR];

  /*! Tangent point latitude [deg]. */
  double tplat[NR];

  /*! Transmittance of ray path. */
  double tau[ND][NR];

  /*! Radiance [W/(m^2 sr cm^-1)]. */
  double rad[ND][NR];

} obs_t;

/**
 * @brief Emissivity look-up tables.
 *
 * Stores precomputed emissivity and source-function data for
 * different gases, spectral channels, and emitter column densities.
 */
typedef struct {

  /*! Number of pressure levels. */
  int np[ND][NG];

  /*! Number of temperatures. */
  int nt[ND][NG][TBLNP];

  /*! Number of column densities. */
  int nu[ND][NG][TBLNP][TBLNT];

  /*! Pressure [hPa]. */
  double p[ND][NG][TBLNP];

  /*! Temperature [K]. */
  double t[ND][NG][TBLNP][TBLNT];

  /*! Column density [molecules/cm^2]. */
  float u[ND][NG][TBLNP][TBLNT][TBLNU];

  /*! Emissivity. */
  float eps[ND][NG][TBLNP][TBLNT][TBLNU];

  /*! Source function temperature [K]. */
  double st[TBLNS];

  /*! Source function radiance [W/(m^2 sr cm^-1)]. */
  double sr[TBLNS][ND];

} tbl_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/**
 * @brief Convert atmospheric data to state vector elements.
 *
 * Extracts selected quantities from an atmospheric profile (@ref atm_t)
 * according to retrieval settings in @ref ctl_t, and appends them to
 * the state vector @p x. For each included quantity, the function also
 * stores its quantity index (@p iqa) and profile index (@p ipa).
 *
 * The function respects retrieval altitude limits defined in @p ctl
 * (e.g., `retp_zmin/zmax`, `rett_zmin/zmax`, etc.) and includes only
 * variables flagged for retrieval (e.g., `ret_clz`, `ret_sft`, etc.).
 *
 * @param[in]  ctl  Control settings defining retrieval configuration and limits.
 * @param[in]  atm  Atmospheric profile data to extract from.
 * @param[out] x    GSL vector to store state-vector elements.
 * @param[out] iqa  Quantity index array corresponding to elements in @p x.
 * @param[out] ipa  Profile index array corresponding to elements in @p x.
 *
 * @return Number of elements written to the state vector.
 *
 * @note Internally calls @ref atm2x_help() to append individual values.
 *
 * @see atm_t, ctl_t, atm2x_help
 *
 * @author Lars Hoffmann
 */
size_t atm2x(
  const ctl_t * ctl,
  const atm_t * atm,
  gsl_vector * x,
  int *iqa,
  int *ipa);

/*! Add element to state vector. */
void atm2x_help(
  const double value,
  const int value_iqa,
  const int value_ip,
  gsl_vector * x,
  int *iqa,
  int *ipa,
  size_t *n);

/*! Convert Cartesian coordinates to geolocation. */
void cart2geo(
  const double *x,
  double *z,
  double *lon,
  double *lat);

/*! Interpolate climatological data. */
void climatology(
  const ctl_t * ctl,
  atm_t * atm_mean);

/*! Compute carbon dioxide continuum (optical depth). */
double ctmco2(
  const double nu,
  const double p,
  const double t,
  const double u);

/*! Compute water vapor continuum (optical depth). */
double ctmh2o(
  const double nu,
  const double p,
  const double t,
  const double q,
  const double u);

/*! Compute nitrogen continuum (absorption coefficient). */
double ctmn2(
  const double nu,
  const double p,
  const double t);

/*! Compute oxygen continuum (absorption coefficient). */
double ctmo2(
  const double nu,
  const double p,
  const double t);

/*! Copy and initialize atmospheric data. */
void copy_atm(
  const ctl_t * ctl,
  atm_t * atm_dest,
  const atm_t * atm_src,
  const int init);

/*! Copy and initialize observation data. */
void copy_obs(
  const ctl_t * ctl,
  obs_t * obs_dest,
  const obs_t * obs_src,
  const int init);

/*! Find index of an emitter. */
int find_emitter(
  const ctl_t * ctl,
  const char *emitter);

/*! Determine ray paths and compute radiative transfer. */
void formod(
  const ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs);

/*! Compute absorption coefficient of continua. */
void formod_continua(
  const ctl_t * ctl,
  const los_t * los,
  const int ip,
  double *beta);

/*! Apply field of view convolution. */
void formod_fov(
  const ctl_t * ctl,
  obs_t * obs);

/*! Compute radiative transfer for a pencil beam. */
void formod_pencil(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const atm_t * atm,
  obs_t * obs,
  const int ir);

/*! Apply RFM for radiative transfer calculations. */
void formod_rfm(
  const ctl_t * ctl,
  const atm_t * atm,
  obs_t * obs);

/*! Compute Planck source function. */
void formod_srcfunc(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const double t,
  double *src);

/*! Convert geolocation to Cartesian coordinates. */
void geo2cart(
  const double z,
  const double lon,
  const double lat,
  double *x);

/*! Set hydrostatic equilibrium. */
void hydrostatic(
  const ctl_t * ctl,
  atm_t * atm);

/*! Determine name of state vector quantity for given index. */
void idx2name(
  const ctl_t * ctl,
  const int idx,
  char *quantity);

/*! Initialize source function table. */
void init_srcfunc(
  const ctl_t * ctl,
  tbl_t * tbl);

/*! Interpolate atmospheric data. */
void intpol_atm(
  const ctl_t * ctl,
  const atm_t * atm,
  const double z,
  double *p,
  double *t,
  double *q,
  double *k);

/*! Get transmittance from look-up tables (CGA method). */
void intpol_tbl_cga(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const los_t * los,
  const int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/*! Get transmittance from look-up tables (EGA method). */
void intpol_tbl_ega(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const los_t * los,
  const int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/*! Interpolate emissivity from look-up tables. */
double intpol_tbl_eps(
  const tbl_t * tbl,
  const int ig,
  const int id,
  const int ip,
  const int it,
  const double u);

/*! Interpolate column density from look-up tables. */
double intpol_tbl_u(
  const tbl_t * tbl,
  const int ig,
  const int id,
  const int ip,
  const int it,
  const double eps);

/*! Convert seconds to date. */
void jsec2time(
  const double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain);

/*! Compute Jacobians. */
void kernel(
  const ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs,
  gsl_matrix * k);

/*! Find array index for irregular grid. */
int locate_irr(
  const double *xx,
  const int n,
  const double x);

/*! Find array index for regular grid. */
int locate_reg(
  const double *xx,
  const int n,
  const double x);

/*! Find array index in float array. */
int locate_tbl(
  const float *xx,
  const int n,
  const double x);

/*! Compose measurement vector. */
size_t obs2y(
  const ctl_t * ctl,
  const obs_t * obs,
  gsl_vector * y,
  int *ida,
  int *ira);

/*! Do ray-tracing to determine LOS. */
void raytrace(
  const ctl_t * ctl,
  const atm_t * atm,
  obs_t * obs,
  los_t * los,
  const int ir);

/*! Read atmospheric data. */
void read_atm(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  atm_t * atm);

/*! Read forward model control parameters. */
void read_ctl(
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read matrix. */
void read_matrix(
  const char *dirname,
  const char *filename,
  gsl_matrix * matrix);

/*! Read observation data. */
void read_obs(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  obs_t * obs);

/*! Read observation data in RFM format. */
double read_obs_rfm(
  const char *basename,
  const double z,
  double *nu,
  double *f,
  int n);

/*! Read RFM spectrum. */
void read_rfm_spec(
  const char *filename,
  double *nu,
  double *rad,
  int *npts);

/*! Read shape function. */
void read_shape(
  const char *filename,
  double *x,
  double *y,
  int *n);

/*! Read look-up table data. */
tbl_t *read_tbl(
  const ctl_t * ctl);

/*! Search control parameter file for variable entry. */
double scan_ctl(
  int argc,
  char *argv[],
  const char *varname,
  const int arridx,
  const char *defvalue,
  char *value);

/*! Calculate solar zenith angle. */
double sza(
  double sec,
  double lon,
  double lat);

/*! Find tangent point of a given LOS. */
void tangent_point(
  const los_t * los,
  double *tpz,
  double *tplon,
  double *tplat);

/*! Convert date to seconds. */
void time2jsec(
  const int year,
  const int mon,
  const int day,
  const int hour,
  const int min,
  const int sec,
  const double remain,
  double *jsec);

/*! Measure wall-clock time. */
void timer(
  const char *name,
  const char *file,
  const char *func,
  int line,
  int mode);

/*! Write atmospheric data. */
void write_atm(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/*! Write atmospheric data in RFM format. */
void write_atm_rfm(
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/*! Write matrix. */
void write_matrix(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const gsl_matrix * matrix,
  const atm_t * atm,
  const obs_t * obs,
  const char *rowspace,
  const char *colspace,
  const char *sort);

/*! Write observation data. */
void write_obs(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const obs_t * obs);

/*! Write shape function. */
void write_shape(
  const char *filename,
  const double *x,
  const double *y,
  const int n);

/*! Write look-up table data. */
void write_tbl(
  const ctl_t * ctl,
  const tbl_t * tbl);

/*! Decompose parameter vector or state vector. */
void x2atm(
  const ctl_t * ctl,
  const gsl_vector * x,
  atm_t * atm);

/*! Get element from state vector. */
void x2atm_help(
  double *value,
  const gsl_vector * x,
  size_t *n);

/*! Decompose measurement vector. */
void y2obs(
  const ctl_t * ctl,
  const gsl_vector * y,
  obs_t * obs);

#endif
