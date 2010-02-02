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

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#if _XOPEN_SOURCE >= 600 || defined (_ISOC99_SOURCE)
#define ISFINITE(x) isfinite (x)
#else
#define ISFINITE(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

#include "cairoint.h"
#include "cairo-spline-offset.h"

static void
_lerp (const point_t *a, const point_t *b, point_t *result, double t)
{
	result->x = (1-t)*a->x + t*b->x;
	result->y = (1-t)*a->y + t*b->y;
}

/* returns false if there could be a discontinuity
 * between the split curves */
static cairo_bool_t
_de_casteljau_t (knots_t *k1, knots_t *k2, double t)
{
	point_t ab, bc, cd;
	point_t abbc, bccd;
	point_t final;

	_lerp (&k1->a, &k1->b, &ab, t);
	_lerp (&k1->b, &k1->c, &bc, t);
	_lerp (&k1->c, &k1->d, &cd, t);
	_lerp (&ab, &bc, &abbc, t);
	_lerp (&bc, &cd, &bccd, t);
	_lerp (&abbc, &bccd, &final, t);

	k2->a = final;
	k2->b = bccd;
	k2->c = cd;
	k2->d = k1->d;

	k1->b = ab;
	k1->c = abbc;
	k1->d = final;

    /* we can test whether the split curve will have a cusp
     * inbetween the two parts by checking the length of the
     * line abbc-bccd. If this line has a non-zero length
     * then the two curves will share a common normal
     * at the points where they meet. Therefore the
     * length must be zero for the normals to be different.
     * However this condition is necessary but not sufficient.
     * XXX: does this matter and how should we deal with it if it does? */
    return abbc.x != final.x || abbc.y != final.y;
}

typedef point_t vector_t;

static vector_t
begin_tangent (knots_t c)
{
    vector_t t;
    if (c.a.x != c.b.x || c.a.y != c.b.y) {
	t.x = c.b.x - c.a.x;
	t.y = c.b.y - c.a.y;
	return t;
    }
    if (c.a.x != c.c.x || c.a.y != c.c.y) {
	t.x = c.c.x - c.a.x;
	t.y = c.c.y - c.a.y;
	return t;
    }
    t.x = c.d.x - c.a.x;
    t.y = c.d.y - c.a.y;
    return t;
}

static vector_t
end_tangent (knots_t c)
{
    vector_t t;
    if (c.d.x != c.c.x || c.d.y != c.c.y) {
	t.x = c.d.x - c.c.x;
	t.y = c.d.y - c.c.y;
	return t;
    }
    if (c.d.x != c.b.x || c.d.y != c.b.y) {
	t.x = c.d.x - c.b.x;
	t.y = c.d.y - c.b.y;
	return t;
    }
    t.x = c.d.x - c.a.x;
    t.y = c.d.y - c.a.y;
    return t;
}

static void
ensure_finite (knots_t k)
{

    if (!(
	    isfinite(k.a.x) &&
	    isfinite(k.a.y) &&
	    isfinite(k.b.x) &&
	    isfinite(k.b.y) &&
	    isfinite(k.c.x) &&
	    isfinite(k.c.y) &&
	    isfinite(k.d.x) &&
	    isfinite(k.d.y))) {
	print_knot(k);
	assert(0);
    }
}

#if 0
static void ensure_smoothness (knots_t a, knots_t b)
{
    double dist_len = hypot(a.d.x - b.a.x, a.d.y - b.a.y);
    print_knot(a);
    print_knot(b);
    printf("dist: %f\n", dist_len);
    vector_t in_norm = end_tangent(a);
    vector_t out_norm = begin_tangent(b);
    double in_a = atan2 (in_norm.y, in_norm.x);
    double out_a = atan2 (out_norm.y, out_norm.x);
    if (in_a < out_a)
	in_a += 2*M_PI;
    double diff = in_a - out_a;
    if (diff > M_PI)
	diff -= M_PI*2;
    if (!(fabs(diff) < 0.001) && !(fabs(diff - M_PI) < 0.001) && !(fabs(diff + M_PI) < 0.001)) {
	/* The angle difference between normals can be either 0 or +-M_PI for a change in direction.
	 * The change in direction can occur in a region of high curvature with a large offset direction.
	 * We move in reverse if the offset distance exceeds the radius of curvature. Note: that you need to
	 * use signed curvature for this too work out properly. Because we only encounter this phenomenom in
	 * the convex areas of the curve.
	 *
	 * XXX: explain why... */
	printf("angle diff: %f %f %f\n", diff, in_a, out_a);
	print_knot(a);
	print_knot(b);
	assert(0);
    }
    if (fabs(dist_len) > 0.01) {
	printf("distlen: %f\n", dist_len);
	print_knot(a);
	print_knot(b);
	assert(0);
    }
}

static void
ensure_smoothness2 (knots_t a, knots_t b, cairo_bool_t ccw) {
    vector_t in_norm = end_tangent(a);
    vector_t out_norm = begin_tangent(b);
    if (!ccw) {
	in_norm = end_tangent(a);
	out_norm = begin_tangent(b);
    }
    double in_a = atan2 (in_norm.y, in_norm.x);
    double out_a = atan2 (out_norm.y, out_norm.x);
    if (in_a < out_a)
	in_a += 2*M_PI;
    double diff = in_a - out_a;
    if (diff > M_PI)
	diff -= M_PI*2;
    if (fabs(diff) > 0.01) {
	printf("gangle diff: %f\n", diff);
	print_knot(a);
	print_knot(b);
    }
}
#endif

static inline vector_t
flip (vector_t a)
{
    vector_t ar;
    ar.x = -a.x;
    ar.y = -a.y;
    return ar;
}

static vector_t
perp (vector_t v)
{
    vector_t p = {-v.y, v.x};
    return p;
}

static inline double
dot (vector_t a, vector_t b)
{
    return a.x * b.x + a.y * b.y;
}

static vector_t
normalize (vector_t v)
{
    double len = hypot(v.x, v.y);
    v.x /= len;
    v.y /= len;
    return v;
}

/* given a normal rotate the vector 90 degrees to the right clockwise
 * This function has a period of 4. e.g. swap(swap(swap(swap(x) == x */
static vector_t
swap (vector_t a)
{
    vector_t ar;
    /* one of these needs to be negative. We choose a.x so that we rotate to the right instead of negating */
    ar.x = a.y;
    ar.y = -a.x;
    return ar;
}

static vector_t
unperp (vector_t a)
{
    return swap (a);
}

typedef struct {
    double t1;
    double t2;
} inflection_points_t;

/* Find the inflection points for a knots_t.
 * The lower inflection point is returned in t1, the higher one in t2.
 *
 * This method for computing inflections points is from:
 * "Fast, precise flattening of cubic Bezier path and offset curves", Hain et. al */
static inflection_points_t
inflection_points (knots_t s)
{
    inflection_points_t ret;

    point_t a = {x:  -s.a.x + 3*s.b.x - 3*s.c.x + s.d.x,  -s.a.y + 3*s.b.y - 3*s.c.y + s.d.y};
    point_t b = {x: 3*s.a.x - 6*s.b.x + 3*s.c.x,         3*s.a.y - 6*s.b.y + 3*s.c.y};
    point_t c = {x:-3*s.a.x + 3*s.b.x,                  -3*s.a.y + 3*s.b.y};
    /* we don't need 'd'
    point_t d = {x:   s.a.x,                               s.a.y};
    */

    /* we're solving for the roots of this equation:
        dx/dt * d2y/dt2 - d2x/dt2 * dy/dt
        == 6*(a.y*b.x - a.x*b.y)*t^2 +6*(a.y*c.x-a.x*c.y)*t + 2*(b.y*c.x - b.x*c.y)
    */
    double t_cusp = (-1./2)*((a.y*c.x - a.x*c.y)/(a.y*b.x - a.x*b.y));
    double t_inner = (t_cusp*t_cusp - ((b.y*c.x - b.x*c.y)/(a.y*b.x - a.x*b.y))/3);

    if (t_inner > 0) {
	double t_adj = sqrt (t_inner);
	ret.t1 = t_cusp - t_adj;
	ret.t2 = t_cusp + t_adj;
    } else {
	ret.t1 = ret.t2 = t_cusp;
    }
    return ret;
}

/* this function could use some tests.
   I believe it only works for angles from 0 to pi. It will always use the smaller arc between the two vectors? */
static knots_t
arc_segment_cart (point_t start, point_t end, double radius, vector_t a, vector_t b)
{
    knots_t arc;
    double r_sin_A = radius * a.y;
    double r_cos_A = radius * a.x;
    double r_sin_B = radius * b.y;
    double r_cos_B = radius * b.x;

    point_t mid, mid2;
    double l, h;

    mid.x = a.x + b.x;
    mid.y = a.y + b.y;
    if (dot (a, b) < 0) {
	/* otherwise, we can flip a, add it
	 * and then use the perpendicular of the result */
	a = flip (a);
	mid.x = a.x + b.x;
	mid.y = a.y + b.y;
	if (dot (a, perp(b)) <= 0) {
	    mid = unperp (mid);
	} else {
	    mid = perp(mid);
	}
    }
    /* normalize */
    l = sqrt(mid.x*mid.x + mid.y*mid.y);
    mid.x /= l;
    mid.y /= l;

    mid2.x = a.x + mid.x;
    mid2.y = a.y + mid.y;

    h = (4./3.)*dot(perp(a), mid2)/dot(a, mid2);

    arc.a.x = start.x;
    arc.a.y = start.y;

    arc.b.x = start.x - h * r_sin_A;
    arc.b.y = start.y + h * r_cos_A;

    arc.c.x =   end.x + h * r_sin_B;
    arc.c.y =   end.y - h * r_cos_B;

    arc.d.x = end.x;
    arc.d.y = end.y;

    return arc;
}

static knots_t
quarter_circle (double xc, double yc, double radius, vector_t normal)
{
    double r_sin_A, r_cos_A;
    double r_sin_B, r_cos_B;
    double h;
    knots_t arc;

    r_sin_A = radius * normal.y;
    r_cos_A = radius * normal.x;
    r_sin_B = radius * -normal.x; /* sin (angle_B); */
    r_cos_B = radius * normal.y;  /* cos (angle_B); */

    h = -0.55228474983079345; /* 4.0/3.0 * tan ((-M_PI/2) / 4.0); */

    arc.a.x = xc + r_cos_A;
    arc.a.y = yc + r_sin_A;

    arc.b.x = xc + r_cos_A - h * r_sin_A;
    arc.b.y = yc + r_sin_A + h * r_cos_A;

    arc.c.x = xc + r_cos_B + h * r_sin_B;
    arc.c.y = yc + r_sin_B - h * r_cos_B;

    arc.d.x = xc + r_cos_B;
    arc.d.y = yc + r_sin_B;

    return arc;
}

static void semi_circle (double xc, double yc, double radius, vector_t out_normal, void (*curve_fn)(void *, knots_t), void *closure)
{
    double len = hypot (out_normal.y, out_normal.x);
    vector_t mid_normal;

    out_normal.x /= len;
    out_normal.y /= len;

    mid_normal.x = out_normal.y; /* we need to rotate by -M_PI/2 */
    mid_normal.y = -out_normal.x;

    curve_fn(closure, quarter_circle (xc, yc, radius, out_normal));
    curve_fn(closure, quarter_circle (xc, yc, radius, mid_normal));
}

static cairo_bool_t
is_knot_flat (knots_t c)
{
    /*
       A well-known flatness test is both cheaper and more reliable
       than the ones you have tried. The essential observation is that
       when the curve is a uniform speed straight line from end to end,
       the control points are evenly spaced from beginning to end.
       Therefore, our measure of how far we deviate from that ideal
       uses distance of the middle controls, not from the line itself,
       but from their ideal *arrangement*. Point 2 should be halfway
       between points 1 and 3; point 3 should be halfway between points
       2 and 4.
       This, too, can be improved. Yes, we can eliminate the square roots
       in the distance tests, retaining the Euclidean metric; but the
       taxicab metric is faster, and also safe. The length of displacement
       (x,y) in this metric is |x|+|y|.
    */

    /* XXX: this isn't necessarily the best test because it
       tests for flatness instead of smallness.
       however, flat things will be offseted easily so we sort of implicitly do the right thing
    */
    /* this value should be much less than the tolerance because we only want to deal with situations where
     * the curvature is extreme. */
    double stop_value = 0.05;
    double error = 0;
    error += fabs((c.c.x - c.b.x) - (c.b.x - c.a.x));
    error += fabs((c.c.y - c.b.y) - (c.b.y - c.a.y));
    error += fabs((c.c.x - c.b.x) - (c.d.x - c.c.x));
    error += fabs((c.c.y - c.b.y) - (c.d.y - c.c.y));
    return error <= stop_value;
}

/* Finds the coefficient with the maximum error. This should be
 * the coefficient closest to area of maximum error.
 * This returns an upper bound of the error. */
static cairo_bool_t
is_curvy_t (double *bezier, int degree, double stop_value, double *t)
{
	double highest_error = .5;
	double error = 0;
	int i=0;
	assert(degree == 6);
	for (i=0; i<=degree; i++) {
	    if (fabs(bezier[i]) > error) {
		error = fabs(bezier[i]);
		highest_error = i/(double)degree;
	    }
	}
	*t = highest_error;

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
static void
approximate_offset_curve_with_shin (knots_t self, double dist, void (*curve_fn)(void *, knots_t), void *closure)
{
    /* since we're going to approximating the offset curve with arcs
       we want to ensure that we don't have any inflection points in our spline */
    inflection_points_t t_inflect = inflection_points(self);

    /* split the curve at the inflection points and convolve each with an approximation of
     * a circle */
    knots_t remaining = self;
    double adjusted_inflect = t_inflect.t2;

    /* check that the first inflection point is in range */
    if (t_inflect.t1 > 0 && t_inflect.t1 < 1) {
	    /* split at the inflection point */
	    _de_casteljau_t (&self, &remaining, t_inflect.t1);

	    /* approximate */
	    curve_fn(closure, convolve_with_circle(self, dist));

	    /* reparamaterize the second inflection point over the 'remaining' spline:
	     * (domain between inflection points)/(domain of remaining) */
	    adjusted_inflect = (t_inflect.t2-t_inflect.t1)/(1 - t_inflect.t1);
    }

    /* check that the second inflection point is in range and not equal to t1.
     * we don't use the reparameterized value so that test remains simple */
    if (t_inflect.t2 > t_inflect.t1 && t_inflect.t2 > 0 && t_inflect.t2 < 1) {
	    /* split into one segment */
	    knots_t self = remaining;
	    _de_casteljau_t (&self, &remaining, adjusted_inflect);
	    curve_fn(closure, convolve_with_circle(self, dist));
    }

    /* deal with what's left of the spline */
    curve_fn(closure, convolve_with_circle(remaining, dist));
}

/* Computes the offset polygon for the control polygon of spline. We can
 * use this is an approximation of the offset spline. This approximation is
 * give by Tiller W, Hanson E G. in "Offsets of two-dimensional profile" */
static knots_t
knot_offset (knots_t self, double width, cairo_bool_t *is_parallel)
{
    /* this code is inspired by the libart mitreing code */

    /* difference vectors between each point in the knot */
    double dx0, dy0, dx1, dy1, dx2, dy2;
    /* difference vector lengths */
    double ld0, ld1, ld2;
    /* temporary vector for dealing with degeneracies */
    double last_x, last_y;

    double scale;

    /* normalized normal vectors for each difference vector */
    double dlx0, dly0,  dlx1, dly1,  dlx2, dly2;

    /* bisected/midpoint vectors */
    double dm1x, dm1y,  dm2x, dm2y;

    double dmr2_1, dmr2_2;

    knots_t result;

    dx0 = self.b.x - self.a.x;
    dy0 = self.b.y - self.a.y;
    ld0 = sqrt(dx0*dx0 + dy0*dy0);

    dx1 = self.c.x - self.b.x;
    dy1 = self.c.y - self.b.y;
    ld1 = sqrt(dx1*dx1 + dy1*dy1);

    dx2 = self.d.x - self.c.x;
    dy2 = self.d.y - self.c.y;
    ld2 = sqrt(dx2*dx2 + dy2*dy2);

    /* compute normalized normals */

    /* deal with any degeneracies in the spline by treating
     * degenerate segments as having the same normal
     * as their neighbour. If all of the segments are degenerate
     * than we will fail the is_parallel test below */
    last_x = 0;
    last_y = 0;

    /* first pass */
    if (ld0) {
	scale = 1. / ld0;
	dlx0 = -dy0 * scale;
	dly0 = dx0 * scale;
	last_x = dlx0;
	last_y = dly0;
    }
    if (ld1) {
	scale = 1. / ld1;
	dlx1 = -dy1 * scale;
	dly1 = dx1 * scale;
	last_x = dlx1;
	last_y = dly1;
    }
    if (ld2) {
	scale = 1. / ld2;
	dlx2 = -dy2 * scale;
	dly2 = dx2 * scale;
	last_x = dlx2;
	last_y = dly2;
    }

    /* second pass */
    if (!ld2) {
	dlx2 = last_x;
	dly2 = last_y;
    } else {
	scale = 1. / ld2;
	dlx2 = -dy2 * scale;
	dly2 = dx2 * scale;
	last_x = dlx2;
	last_y = dly2;
    }
    if (!ld1) {
	dlx1 = last_x;
	dly1 = last_y;
    } else {
	scale = 1. / ld1;
	dlx1 = -dy1 * scale;
	dly1 = dx1 * scale;
	last_x = dlx1;
	last_y = dly1;
    }
    if (!ld0) {
	dlx0 = last_x;
	dly0 = last_y;
    } else {
	scale = 1. / ld0;
	dlx0 = -dy0 * scale;
	dly0 = dx0 * scale;
	last_x = dlx0;
	last_y = dly0;
    }

    /* mid-point vector sum */
    dm1x = (dlx0 + dlx1);
    dm1y = (dly0 + dly1);

    dm2x = (dlx1 + dlx2);
    dm2y = (dly1 + dly2);

    /* length squared of the mid-point vector sum */
    dmr2_1 = (dm1x * dm1x + dm1y * dm1y);
    dmr2_2 = (dm2x * dm2x + dm2y * dm2y);

    /* the choice of this EPSILON is arbitrary and could use further study */
#define PARALLEL_EPSILON 1e-6
    if (fabs(dmr2_1) < PARALLEL_EPSILON || fabs(dmr2_2) < PARALLEL_EPSILON)
	*is_parallel = TRUE;
    else
	*is_parallel = FALSE;

    scale = width * 2 / dmr2_1;
    dm1x *= scale;
    dm1y *= scale;

    scale = width * 2 / dmr2_2;
    dm2x *= scale;
    dm2y *= scale;

    /* write out the result */
    result.a.x = self.a.x + dlx0 * width;
    result.a.y = self.a.y + dly0 * width;

    result.b.x = self.b.x + dm1x;
    result.b.y = self.b.y + dm1y;

    result.c.x = self.c.x + dm2x;
    result.c.y = self.c.y + dm2y;

    result.d.x = self.d.x + dlx2 * width;
    result.d.y = self.d.y + dly2 * width;

    return result;
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
void
curve_offset (knots_t self, double dist, double tolerance, void (*curve_fn)(void *, knots_t), void *closure)
{
    cairo_bool_t recurse;
    double break_point;
    double max_error = tolerance;
    knots_t offset = knot_offset(self, dist, &recurse);

    if (!recurse) {
	/* we need to make sure we have a finite offset before checking the error */
	ensure_finite(offset);
	recurse = !error_within_bounds_elber(self, offset, dist, max_error, &break_point);
    } else {
	/* TODO: we could probably choose a better place to split than halfway
	 * for times when we know we're going to recurse */
	break_point = .5;
    }

    if (recurse) {
	/* error is too large. subdivide on max t and try again. */
	knots_t left = self;
	knots_t right;

	cairo_bool_t is_continuous;

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
	if (is_knot_flat(self)) {
	    /* we don't want to recurse anymore because we're pretty flat,
	     * instead approximate the offset curve with an arc.
	     *
	     * The previously generated offset curve can become really bad if the the intersections
	     * are far away from the original curve (self) */

	    approximate_offset_curve_with_shin (self, dist, curve_fn, closure);
	    return;
	}

	is_continuous = _de_casteljau_t(&left, &right, break_point);

	curve_offset(left, dist, tolerance, curve_fn, closure);

	if (!is_continuous) {
	    /* add circle:
	       we can use the exit the normal from the left side
	       as the begining normal of our cusp */
	    vector_t normal = perp (end_tangent (left));
	    semi_circle (left.d.x, left.d.y, dist, normal, curve_fn, closure);
	}

	curve_offset(right, dist, tolerance, curve_fn, closure);
    } else {
	curve_fn(closure, offset);
    }
}

