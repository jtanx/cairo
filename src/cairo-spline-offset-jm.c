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

#define _GNU_SOURCE

#include "cairoint.h"

#include "cairo-spline-offset.h"

/* Finds the coefficient with the maximum error. This should be
 * the coefficient closest to area of maximum error.
 * This returns an upper bound of the error. */
static double
error_upper_bound (const double bezier[7], int degree, double stop_value, double *t)
{
    double worst_t = .5;
    double max_error = 0;
    int i;

    for (i=0; i<7; i++)
	if (fabs(bezier[i]) > error) {
	    error = fabs(bezier[i]);
	    worst_t = i/6.;
	}

    *t = worst_t;

	/* we shouldn't ever split at the edges because the approximation
	 * should be perfect there. XXX: we might be able to take
	 * advantage of this when computing the error polynomial. */
#if 0
	if (*t == 0 || *t == 1) {
	    for (i=0; i<=degree; i++) {
		error = fabs(bezier[i]);
		printf("error: %f\n", error);
	    }
	    assert(0);
	}
#endif
	/* *t = 0.5; */
	return error >= stop_value;
}

static cairo_bool_t
error_within_bounds_elber (knots_t bz1, knots_t bz_offset_approx,
		double desired_offset,
		double max_error,
		double *t)
{
    double bzprod[3+3 +1] = {};

    /* Elber and Cohnen propose the following method of computing the
       error between an curve and it's offset:
       e(t) = || C(t) - Coffset_approx(t) ||^2 - offset^2

       It is important to note that this function doesn't represent a unique offset curve.
       For example, it's possible for the offset curve to cross the original spline
       provided it stays the appropriate distance away:

                 xxx
                x
               ****
             **x
            *  x
           *   x
            *  x
             **
               ****
               x
                xxx

      "Planar Curve Offset Base on Circle Approximation" suggests further problems
      with this error measure. They give the following example: C(t) + (r, 0) will
      have an error mesaure of 0, but is not even close to the actual offset curve.

      This paper suggests using the angle between the difference vector and the normal
      at point t, but doesn't go into much detail.
     */

    /* Now we compute this error function.
     * First, subtract the approx */
    bz1.a.x -= bz_offset_approx.a.x;
    bz1.a.y -= bz_offset_approx.a.y;
    bz1.b.x -= bz_offset_approx.b.x;
    bz1.b.y -= bz_offset_approx.b.y;
    bz1.c.x -= bz_offset_approx.c.x;
    bz1.c.y -= bz_offset_approx.c.y;
    bz1.d.x -= bz_offset_approx.d.x;
    bz1.d.y -= bz_offset_approx.d.y;
/*
     Next we compute the dot product of bz1 with itself This involves squaring
     the x and y polynomials and adding the result together.

     However, we want to multiply the coefficients so that we can treat the
     result as a higher degree bezier curve.

     The following method is given in the Symbolic Computation section of
     "Comparing Offset Curve Approximation Methods" by Elber, Lee and Kim

  k     l       ⎛ k⎞  i     k-i⎛ l⎞  j       l-j
 B (t) B (t) =  ⎜  ⎟ t (1-t)   ⎜  ⎟ t   (1-t)
  i     j       ⎝ i⎠           ⎝ j⎠
	        ⎛ k⎞ ⎛ l⎞  i+j     k+l-i-j
	     =  ⎜  ⎟ ⎜  ⎟ t   (1-t)
	        ⎝ i⎠ ⎝ j⎠

	        ⎛ k⎞ ⎛ l⎞
	     =  ⎝ i⎠ ⎝ j⎠  k+l
	         ──────── B   (t)
	        ⎛ k + l⎞   i+j
	        ⎝ i + j⎠

               m               n
C1(t) C2(t) =  ∑ P_i B_i,m (t) ∑ P_j B_j,n (t)
              i=0             j=0

	       m   n
            =  ∑   ∑ P_i P_j B_i,m (t) B_j,n (t)
              i=0 j=0

	       m+n       m + n
	    =   ∑ Q_k B_k     (t)
	       k=0

Where Q_k is:
		       ⎛ m⎞ ⎛ n⎞
		       ⎝ i⎠ ⎝ j⎠
       Q_k = P_i * P_j ─────────
			⎛ m+n⎞
			⎝ i+j⎠

And for the case of degree 3 beziers

		       ⎛ 3⎞ ⎛ 3⎞
		       ⎝ i⎠ ⎝ j⎠
       Q_k = P_i * P_j ─────────
			⎛  6 ⎞
			⎝ i+j⎠

      This expands to:

	i j
	0 0 0 = (1 * 1) / 1  = 1

	0 1 1 = (1 * 3) / 6  = .5
	1 0 1 = (3 * 1) / 6  = .5

	0 2 2 = (1 * 3) / 15 = .2
	1 1 2 = (3 * 3) / 15 = .6
	2 0 2 = (3 * 1) / 15 = .2

	0 3 3 = (1 * 1) / 20 = .05
	1 2 3 = (3 * 3) / 20 = .45
	2 1 3 = (3 * 3) / 20 = .45
	3 0 3 = (1 * 1) / 20 = .05

	1 3 4 = (3 * 1) / 15 = .2
	2 2 4 = (3 * 3) / 15 = .6
	3 1 4 = (1 * 3) / 15 = .2

	2 3 5 = (3 * 1) / 6  = .5
	3 2 5 = (1 * 3) / 6  = .5

	3 3 6 = (1 * 1) / 1  = 1

Which simplifies to the following:
	*/


	bzprod[0] += bz1.a.x * bz1.a.x;
	bzprod[0] += bz1.a.y * bz1.a.y;

	bzprod[1] += bz1.a.x * bz1.b.x;
	bzprod[1] += bz1.a.y * bz1.b.y;

	bzprod[2] += bz1.a.x * bz1.c.x * 0.4;
	bzprod[2] += bz1.b.x * bz1.b.x * 0.6;
	bzprod[2] += bz1.a.y * bz1.c.y * 0.4;
	bzprod[2] += bz1.b.y * bz1.b.y * 0.6;

	bzprod[3] += bz1.d.x * bz1.a.x * 0.1;
	bzprod[3] += bz1.b.x * bz1.c.x * 0.9;
	bzprod[3] += bz1.d.y * bz1.a.y * 0.1;
	bzprod[3] += bz1.b.y * bz1.c.y * 0.9;

	bzprod[4] += bz1.b.x * bz1.d.x * 0.4;
	bzprod[4] += bz1.c.x * bz1.c.x * 0.6;
	bzprod[4] += bz1.d.y * bz1.b.y * 0.4;
	bzprod[4] += bz1.c.y * bz1.c.y * 0.6;

	bzprod[5] += bz1.c.x * bz1.d.x;
	bzprod[5] += bz1.d.y * bz1.c.y;

	bzprod[6] += bz1.d.x * bz1.d.x;
	bzprod[6] += bz1.d.y * bz1.d.y;


	/* the result is a bezier polynomial that represents the distance squared between
	 * the two input polynomials */

	/* we don't need to subtract the desired offset because our metric doesn't depend on being centered around 0 */
#if 1
	bzprod[0] -= (desired_offset*desired_offset);
	bzprod[1] -= (desired_offset*desired_offset);
	bzprod[2] -= (desired_offset*desired_offset);
	bzprod[3] -= (desired_offset*desired_offset);
	bzprod[4] -= (desired_offset*desired_offset);
	bzprod[5] -= (desired_offset*desired_offset);
	bzprod[6] -= (desired_offset*desired_offset);
#endif
	*t = 0.5;
	return !is_curvy_t(bzprod, sizeof(bzprod)/sizeof(bzprod[0])-1, max_error, t);
}

void
print_knot (knots_t c) {
	/*printf("(%.13a %.13a) (%.13a %.13a) (%.13a %.13a) (%.13a %.13a)\n", c.a.x, c.a.y, c.b.x, c.b.y, c.c.x, c.c.y, c.d.x, c.d.y);*/
	printf("(%.5f %.5f) (%.5f %.5f) (%.5f %.5f) (%.5f %.5f)\n", c.a.x, c.a.y, c.b.x, c.b.y, c.c.x, c.c.y, c.d.x, c.d.y);
}

static knots_t
convolve_with_circle (knots_t spline, double dist)
{
    vector_t in_normal  = normalize(perp(begin_tangent(spline)));
    vector_t out_normal = normalize(perp(end_tangent(spline)));
    /* find an arc the goes from the input normal to the output_normal */
    knots_t arc_spline = arc_segment_cart(in_normal, out_normal, 1, in_normal, out_normal);

    /* approximate convolving 'spline' with the arc spline (XXX: is convolve the right word?) */
    spline.a.x += dist * arc_spline.a.x;
    spline.a.y += dist * arc_spline.a.y;
    spline.b.x += dist * arc_spline.b.x;
    spline.b.y += dist * arc_spline.b.y;
    spline.c.x += dist * arc_spline.c.x;
    spline.c.y += dist * arc_spline.c.y;
    spline.d.x += dist * arc_spline.d.x;
    spline.d.y += dist * arc_spline.d.y;

    return spline;
}

/* This approach is inspired by "Approximation of circular arcs and offset curves by Bezier curves of high degree"
 * by YJ Ahn, Y soo Kim, Y Shin.
 *
 * This is a similar idea to:
 *  "Offset approximation based on reparameterizaing the path of a moving point along the base circle", ZHAO, WANG
 *  and a bunch of papers by Lee et al.
 * "Planar curve oﬀset based on circle approximation"
 * "New approximation methods of planar oﬀset and convolution curves"
 * "Polynomial/rational approximation of Minkowski sum boundary curves"
 *
 * However these other approaches produce curves of higher degree than the input curve making them unsuitable
 * for use here.
 * */
static cairo_status_t
approximate_offset_curve_with_shin (const cairo_spline_knots_t *self, double dist, cairo_status_t (*curve_fn)(void *, const cairo_spline_knots_t *), void *closure)
{
    /* since we're going to approximating the offset curve with arcs
       we want to ensure that we don't have any inflection points in our spline */
    double t_inflect[2] = _cairo_spline_inflection_points (k);

    /* split the curve at the inflection points and convolve each with an approximation of
     * a circle */
    cairo_spline_knots_t remaining = *self;
    double adjusted_inflect = t_inflect[1];

    /* check that the first inflection point is in range and not equal to t2
     * (we avoid reparametrizing the curve and just do the second split if
     * t1==t2) */
    if (t_inflect[0] > 0 && t_inflect[0] < 1 && t_inflect[0] < t_inflect[1]) {
	/* split at the inflection point */
	cairo_spline_knots_t first;
	_cairo_spline_de_casteljau (&remaining, &first, &remaining, t_inflect[0]);

	/* approximate */
	if (!_cairo_spline_isdegenerate (&first)) {
	    convolve_with_circle (&first, &first, dist);
	    status = curve_fn(closure, &first);
	    if (unlikely (status))
		return status;
	}

	/* reparamaterize the second inflection point over the 'remaining' spline:
	 * (domain between inflection points)/(domain of remaining) */
	t_inflect[1] = (t_inflect[1] - t_inflect[0]) / (1 - t_inflect[0]);
    }

    if (t_inflect[1] > 0 && t_inflect[1] < 1) {
	/* split at the inflection point */
	cairo_spline_knots_t second;
	_cairo_spline_de_casteljau (&remaining, &second, &remaining, t_inflect[1]);

	/* approximate */
	if (!_cairo_spline_isdegenerate (&second)) {
	    convolve_with_circle (&second, &second, dist);
	    status = curve_fn(closure, &second);
	    if (unlikely (status))
		return status;
	}
    }

    /* approximate the last part of the spline */
    if (!_cairo_spline_isdegenerate (&remaining)) {
	convolve_with_circle (&remaining, &remaining, dist);
	status = curve_fn(closure, &remaining);
    }
    return status;
}

/* Computes the offset polygon for the control polygon of spline. We can
 * use this is an approximation of the offset spline. This approximation is
 * give by Tiller W, Hanson E G. in "Offsets of two-dimensional profile" */
static cairo_bool_t
knot_offset (cairo_spline_knots_t *dest, const cairo_spline_knots_t *src, double width)
{
    /* this code is inspired by the libart mitreing code */

    /* difference vectors between each point in the knot */
    cairo_vector_t diff[3];
    /* offset and average offset vectors */
    cairo_vector_t off[3], mid[2];

    diff[0] = _cairo_vector_sub (self->b, self->a);
    diff[1] = _cairo_vector_sub (self->c, self->b);
    diff[2] = _cairo_vector_sub (self->d, self->c);

    /* compute normalized normals */
    /* deal with any degeneracies in the spline by treating
     * degenerate segments as having the same normal
     * as their neighbour. If all of the segments are degenerate
     * than we will fail the is_parallel test below */
    for (i=0; i<3; i++)
	if (!_cairo_vector_iszero (diff[i]))
	    off[i] = _cairo_vector_perpendicular (_cairo_vector_normalize (diff[i]));

    if (_cairo_vector_iszero (diff[1]))
	off[1] = _cairo_vector_iszero(diff[2]) ? off[0] : off[2]

    if (_cairo_vector_iszero (diff[2]))
	off[2] = off[1];

    if (_cairo_vector_iszero (diff[0]))
	off[0] = off[1];

    /* mid-point vector sum */
    mid[0] = _cairo_vector_add (off[0], off[1]);
    mid[1] = _cairo_vector_add (off[1], off[2]);

    /* TODO: check nongeneracy */
    /* the choice of this EPSILON is arbitrary and could use further study */
#define PARALLEL_EPSILON 1e-6
    if (fabs (_cairo_vector_dot(mid[0], mid[0])) < PARALLEL_EPSILON || fabs (_cairo_vector_dot(mid[1], mid[1])) < PARALLEL_EPSILON)
	return FALSE;

    mid[0] = _cairo_vector_normalize (mid[0]);
    mid[1] = _cairo_vector_normalize (mid[1]);

    /* write out the result */
    dest->a = _cairo_vector_add (src->a, _cairo_vector_scale (off[0], width));
    dest->b = _cairo_vector_add (src->b, _cairo_vector_scale (mid[0], 2 * width));
    dest->c = _cairo_vector_add (src->c, _cairo_vector_scale (mid[1], 2 * width));
    dest->d = _cairo_vector_add (src->d, _cairo_vector_scale (off[2], width));

    return TRUE;
}

/* plan:
 * if we approximate the original bezier curve with line segments
 * and then produce an offset to those lines. We will end up
 * with sharp angles at places of high curvature.
 *
 * It is at these places of high curvature that we want to approximate
 * the offset with an arc instead of continuing to subdivide.
 *
 * How do we find the areas of high curvature?
 * - Ideas:
 *   As we subdivide the original bezier
 *   we can get an idea of how flat it is.
 *   If the original bezier is flat but the offset approximation is not within our bounds
 *   we should stop recursing and approximate the segment with an arc.
 *
 *   The degree of required flatness of the original curve should be related to the
 *   the offset length. The greater the offset distance the less flat the original
 *   bezier needs to be.
 *
 *   pros:
 *   - this gives us an easier condition to reason about the termination of the
 *     algorithm because we terminate the recursion when the input bezier has
 *     been subdivided to 'flat' instead of when the offset is good.
 *     This is valuable because it's possible to create an input curve that does not have a
 *     good offset approximation down to the resolution of a double.
 */
cairo_status_t
curve_offset (const cairo_spline_knots_t *k, double dist, double tolerance, cairo_status_t (*curve_fn)(void *, const cairo_spline_knots_t *), void *closure)
{
    cairo_knots_t offset;
    cairo_bool_t recurse = knot_offset(&offset, k, dist);
    double break_point = 0.5;
    /* TODO: we could probably choose a better place to split than halfway
     * for times when we know we're going to recurse */

    if (!recurse) {
	assert (_cairo_spline_isfinite (offset));
	/* error is too large. subdivide on max t and try again. */
	recurse = error_upper_bound (k, &offset, dist, &break_point) > tolerance;
    }

    if (!recurse) {
	return curve_fn (closure, &offset);
    } else  if (_cairo_spline_curvyness (k) < tolerance * 0.125) {
	/* We need to do something to handle regions of high curvature
	 *
	 * skia uses the first and second derivitives at the points of 'max curvature' to check for
	 * pinchynes.
	 * qt tests whether the points are close and reverse direction.
	 *   float l = (orig->x1 - orig->x2)*(orig->x1 - orig->x2) + (orig->y3 - orig->y4)*(orig->y3 - orig->y4);
	 * float dot = (orig->x1 - orig->x2)*(orig->x3 - orig->x4) + (orig->y1 - orig->y2)*(orig->y3 - orig->y4);
	 *
	 * we just check if the knot is flat and has large error.
	 */
	/* we don't want to recurse anymore because we're pretty flat,
	 * instead approximate the offset curve with an arc.
	 *
	 * The previously generated offset curve can become really bad if the the intersections
	 * are far away from the original curve (self) */
	return approximate_offset_curve_with_shin (self, dist, curve_fn, closure);
    } else {
	cairo_bool_t is_continuous = _cairo_spline_de_casteljau (k, &left, &right, break_point);

	status = curve_offset(&left, dist, tolerance, curve_fn, closure);
	if (unlikely (status))
	    return status;

	if (!is_continuous) {
	    /* add circle:
	       we can use the exit the normal from the left side
	       as the begining normal of our cusp */
	    cairo_vector_t mid = _cairo_vector_normalize (_cairo_spline_end_tangent (&left));
	    cairo_vector_t start = _cairo_vector_perpendicular (mid);
	    cairo_vector_t end = _cairo_vector_opposite (start);
	    cairo_spline_knots_t arc;
	    
	    _cairo_spline_arc (&arc, left.d, dist, &start, &mid);
	    status = curve_fn (&arc, closure);
	    if (unlikely (status))
		return status;
	    
	    _cairo_spline_arc (&arc, left.d, dist, &mid, &end);
	    status = curve_fn (&arc, closure);
	    if (unlikely (status))
		return status;
	}

	return curve_offset (&right, dist, tolerance, curve_fn, closure);
    }
}
