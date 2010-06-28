/* -*- Mode: c; tab-width: 8; c-basic-offset: 4; indent-tabs-mode: t; -*- */
/* cairo - a vector graphics library with display and print output
 *
 * Copyright © 2009 Mozilla Foundation
 * Copyright © 2010 Jeff Muizelaar
 *
 * This library is free software; you can redistribute it and/or
 * modify it either under the terms of the GNU Lesser General Public
 * License version 2.1 as published by the Free Software Foundation
 * (the "LGPL") or, at your option, under the terms of the Mozilla
 * Public License Version 1.1 (the "MPL"). If you do not alter this
 * notice, a recipient may use your version of this file under either
 * the MPL or the LGPL.
 *
 * You should have received a copy of the LGPL along with this library
 * in the file COPYING-LGPL-2.1; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * You should have received a copy of the MPL along with this library
 * in the file COPYING-MPL-1.1
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY
 * OF ANY KIND, either express or implied. See the LGPL or the MPL for
 * the specific language governing rights and limitations.
 *
 * The Original Code is the cairo graphics library.
 *
 * The Initial Developer of the Original Code is Jeff Muizelaar
 *
 */

#include "cairo-spline-private.h"
#include "cairo-spline-offset.h"

/* assumes that the input spline is not degenerate and that its normal sweeps less than pi.
   src and dest MUST be different splines */
static void
_cairo_spline_offset_circle (cairo_spline_knots_double_t *dest, const cairo_spline_knots_double_t *src, double dist)
{
    cairo_vector_t in_normal = _cairo_vector_perpendicular (_cairo_spline_begin_tangent (src));
    cairo_vector_t out_normal = _cairo_vector_perpendicular (_cairo_spline_begin_tangent (src));

    assert (src != dest);

    /* find an arc the goes from the input normal to the output_normal */
    _cairo_spline_arc (dest, &_cairo_vector_zero, dist, &in_normal, &out_normal);

    /* approximate convolving 'spline' with the arc spline */
    _cairo_spline_scale (dest, dest, dist);
    _cairo_spline_add (dest, dest, src);
}

static cairo_status_t
_cairo_spline_do_arc (const cairo_vector_t *center, double radius, const cairo_vector_t *start, const cairo_vector_t *end, double tolerance, cairo_status_t (*curve_fn)(void *, const cairo_spline_knots_double_t *), void *closure)
{
    cairo_spline_knots_double_t arc;
    cairo_vector_t tmp, tmp2;
    double maxsqradius;
    int subdivision;
    cairo_status_t status;

    maxsqradius = radius + tolerance;
    maxsqradius *= maxsqradius;

    for (subdivision = 1; _cairo_spline_arc (&arc, center, radius, start, &tmp) > maxsqradius; subdivision <<= 1)
	tmp = _cairo_vector_bisect (*start, tmp);

    switch (subdivision) {
    default:
	for (tmp2 = *start; subdivision > 2; subdivision -= 2) {
	    status = curve_fn (closure, &arc);
	    if (unlikely (status))
		return status;

	    tmp2 = _cairo_vector_add (tmp, _cairo_vector_opposite (_cairo_vector_sub (tmp2, tmp)));
	    _cairo_spline_arc (&arc, center, radius, &tmp, &tmp2);
	    status = curve_fn (closure, &arc);
	    if (unlikely (status))
		return status;

	    tmp = _cairo_vector_add (tmp2, _cairo_vector_opposite (_cairo_vector_sub (tmp, tmp2)));
	    _cairo_spline_arc (&arc, center, radius, &tmp2, &tmp);
	}
	/* fall through */
    case 2:
	status = curve_fn (closure, &arc);
	if (unlikely (status))
	    return status;
	_cairo_spline_arc (&arc, center, radius, &tmp, end);
	/* fall through */
    case 1:
	return curve_fn (closure, &arc);
    }
}

static cairo_bool_t
_cairo_spline_offset_needs_recursion (const cairo_spline_knots_double_t *k, double dist, double tolerance, double *t)
{
   double sweep[2];
 
   _cairo_spline_normal_sweep (k, sweep);

   *t = 0.5;

   return fabs (sweep[0]) > 1.0 || fabs (sweep[1]) > 1.0; /* TODO: additional checks for tolerance */
}

/* k is modified by the function */
static cairo_status_t
_cairo_spline_offset_simple (cairo_spline_knots_double_t *k, double dist, double tolerance, cairo_status_t (*curve_fn)(void *, const cairo_spline_knots_double_t *), void *closure)
{
    cairo_spline_knots_double_t tmp;
    cairo_status_t status;
    double t;

    while (1) {
	if (_cairo_spline_isdegenerate (k))
	    return CAIRO_STATUS_SUCCESS;

	if (_cairo_spline_curvyness (k) < tolerance) {
	    /* TODO: respect tolerance */
	    _cairo_spline_offset_circle (&tmp, k, dist);
	    return curve_fn (closure, &tmp);
	}

	if (_cairo_spline_offset_needs_recursion (k, dist, tolerance, &t)) {
	    _cairo_spline_de_casteljau (k, &tmp, k, t);

	    status = _cairo_spline_offset_simple (&tmp, dist, tolerance, curve_fn, closure);
	    if (unlikely (status))
		return status;

	    /* tail-recursion*/
	} else {
	    _cairo_spline_offset_circle (&tmp, k, dist);
	    return curve_fn (closure, &tmp);
	}
    }
}

cairo_status_t
_cairo_spline_offset (const cairo_spline_knots_double_t *k, double dist, double tolerance, cairo_status_t (*curve_fn)(void *, const cairo_spline_knots_double_t *), void *closure)
{
    cairo_spline_knots_double_t split[3];
    cairo_bool_t cusp[2];
    int i, parts = _cairo_spline_split_inflection (split, cusp, k);
    cairo_status_t status;

    for (i = 0; i < parts - 1; i++) {
	status = _cairo_spline_offset_simple (split+i, dist, tolerance, curve_fn, closure);
	if (unlikely (status))
	    return status;
	if (cusp[i]) {
	    cairo_vector_t begin = _cairo_vector_perpendicular (_cairo_spline_end_tangent (split+i));
	    cairo_vector_t end = _cairo_vector_perpendicular (_cairo_spline_begin_tangent (split+i+1));
	    status = _cairo_spline_do_arc (&split[i].d, dist, &begin, &end, tolerance, curve_fn, closure);
	    if (unlikely (status))
		return status;
	}	    
    }

    return _cairo_spline_offset_simple (split+parts-1, dist, tolerance, curve_fn, closure);;
}
