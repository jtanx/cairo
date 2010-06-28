/* -*- Mode: c; tab-width: 8; c-basic-offset: 4; indent-tabs-mode: t; -*- */
#ifndef CAIRO_SPLINE_PRIVATE_H
#define CAIRO_SPLINE_PRIVATE_H

#include "cairo-vector-private.h"

typedef struct _cairo_spline_knots_double {
    cairo_point_double_t a, b, c, d;
} cairo_spline_knots_double_t;

static cairo_bool_t
_cairo_spline_isdegenerate (const cairo_spline_knots_double_t *k)
{
    return
	_cairo_vector_equals (k->a, k->b) &&
	_cairo_vector_equals (k->a, k->c) &&
	_cairo_vector_equals (k->a, k->d);
}

static cairo_vector_t
_cairo_spline_begin_tangent (const cairo_spline_knots_double_t *c)
{
    cairo_vector_t r;
    assert (!_cairo_spline_isdegenerate (c));

    if (!_cairo_vector_equals (c->a, c->b))
	r = _cairo_vector_sub (c->b, c->a);
    else if (!_cairo_vector_equals (c->a, c->c))
	r = _cairo_vector_sub (c->c, c->a);
    else
	r = _cairo_vector_sub (c->d, c->a);

    return _cairo_vector_normalize (r);
}

static cairo_vector_t
_cairo_spline_end_tangent (const cairo_spline_knots_double_t *c)
{
    cairo_vector_t r;
    assert (!_cairo_spline_isdegenerate (c));

    if (!_cairo_vector_equals (c->d, c->c))
	r = _cairo_vector_sub (c->d, c->c);
    else if (!_cairo_vector_equals (c->d, c->b))
	r = _cairo_vector_sub (c->d, c->b);
    else
	r = _cairo_vector_sub (c->d, c->a);
    return _cairo_vector_normalize (r);
}

static cairo_bool_t
_cairo_spline_isfinite (const cairo_spline_knots_double_t *k)
{
    return
	_cairo_vector_isfinite (k->a) &&
	_cairo_vector_isfinite (k->b) &&
	_cairo_vector_isfinite (k->c) &&
	_cairo_vector_isfinite (k->d);
}

static double
_cairo_spline_curvyness (const cairo_spline_knots_double_t *k)
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
    cairo_vector_t ba = _cairo_vector_sub (k->b, k->a);
    cairo_vector_t cb = _cairo_vector_sub (k->c, k->b);
    cairo_vector_t dc = _cairo_vector_sub (k->d, k->c);
    cairo_vector_t e1 = _cairo_vector_sub (cb, ba);
    cairo_vector_t e2 = _cairo_vector_sub (cb, dc);
    return fabs (e1.x) + fabs (e1.y) + fabs (e2.x) + fabs (e2.y);
}


static cairo_bool_t
_cairo_spline_de_casteljau (const cairo_spline_knots_double_t *k, cairo_spline_knots_double_t *k1, cairo_spline_knots_double_t *k2, double t)
{
    cairo_vector_t tmp = _cairo_vector_lerp (k->b, k->c, t); /* bc */
	
    k1->a = k->a;
    k2->d = k->d;

    k1->b = _cairo_vector_lerp (k->a, k->b, t); /* ab */
    k2->c = _cairo_vector_lerp (k->c, k->d, t); /* cd */

    k1->c = _cairo_vector_lerp (k1->b, tmp, t); /* abbc */
    k2->b = _cairo_vector_lerp (tmp, k2->c, t); /* bccd */

    k1->d = _cairo_vector_lerp (k1->c, k2->b, t);
    k2->a = k1->d;

    /* we can test whether the split curve will have a cusp
     * inbetween the two parts by checking abbc and bccd.
     * A cusp can be present iff they are equal (unless the
     * spline is degenerate) */
    return _cairo_vector_equals (k1->c, k2->b);
}

/* Find the inflection points for a knots_t.
 * The lower inflection point is returned in t1, the higher one in t2.
 *
 * This method for computing inflections points is from:
 * "Fast, precise flattening of cubic Bezier path and offset curves", Hain et. al */
static int
_cairo_spline_inflection_points (const cairo_spline_knots_double_t *k, double *t_inflection)
{
    double den;
    int r = 0;

    cairo_vector_t a = {  -k->a.x + 3*k->b.x - 3*k->c.x + k->d.x,  -k->a.y + 3*k->b.y - 3*k->c.y + k->d.y};
    cairo_vector_t b = { 3*k->a.x - 6*k->b.x + 3*k->c.x,          3*k->a.y - 6*k->b.y + 3*k->c.y};
    cairo_vector_t c = {-3*k->a.x + 3*k->b.x,                    -3*k->a.y + 3*k->b.y};
    /* we don't need 'd'
    cairo_vector_t d = {   k->a.x,                                  k->a.y};
    */

    /* we're solving for the roots of this equation:
        dx/dt * d2y/dt2 - d2x/dt2 * dy/dt
        == 6*(a.y*b.x - a.x*b.y)*t^2 +6*(a.y*c.x-a.x*c.y)*t + 2*(b.y*c.x - b.x*c.y)
    */
    if ((den = _cairo_vector_cross(a, b)) != 0.0) {
	double invden = 1 / den;
	double t_cusp = -0.5 * _cairo_vector_cross (a, c) * invden;
	double t_inner = t_cusp*t_cusp - _cairo_vector_cross (b, c) * invden / 3;
	
	if (t_inner > 0.0) {
	    double t_adj = sqrt (t_inner);
	    t_inflection[r] = t_cusp - t_adj;
	    if (t_inflection[r] > 0.0 && t_inflection[r] < 1.0)
		r++;
	    t_inflection[r] = t_cusp + t_adj;
	} else {
	    t_inflection[r] = t_cusp;
	}
    } else if ((den = _cairo_vector_cross(a, c)) != 0.0) {
	t_inflection[r] = - _cairo_vector_cross (b, c) / (den * 3);
    } else
	return 0;

    if (t_inflection[r] > 0.0 && t_inflection[r] < 1.0)
	r++;
    
    return r;
}

static void
_cairo_spline_translate (cairo_spline_knots_double_t *dest, const cairo_spline_knots_double_t *src, const cairo_vector_t *t)
{
    dest->a = _cairo_vector_add (src->a, *t);
    dest->b = _cairo_vector_add (src->b, *t);
    dest->c = _cairo_vector_add (src->c, *t);
    dest->d = _cairo_vector_add (src->d, *t);
}

static void
_cairo_spline_scale (cairo_spline_knots_double_t *dest, const cairo_spline_knots_double_t *src, double scale)
{
    dest->a = _cairo_vector_scale (src->a, scale);
    dest->b = _cairo_vector_scale (src->b, scale);
    dest->c = _cairo_vector_scale (src->c, scale);
    dest->d = _cairo_vector_scale (src->d, scale);
}

static void
_cairo_spline_add (cairo_spline_knots_double_t *dest, const cairo_spline_knots_double_t *a, const cairo_spline_knots_double_t *b)
{
    dest->a = _cairo_vector_add (a->a, b->a);
    dest->b = _cairo_vector_add (a->b, b->b);
    dest->c = _cairo_vector_add (a->c, b->c);
    dest->d = _cairo_vector_add (a->d, b->d);
}

/* Compute a spline approximation of the arc
   centered at xc, yc from the angle a to the angle b

   The angle between a and b should not be more than a
   quarter circle (pi/2)

   The approximation is similar to an approximation given in:
   "Approximation of a cubic bezier curve by circular arcs and vice versa"
   by Alekas Riškus. However that approximation becomes unstable when the
   angle of the arc approaches 0.

   This approximation is inspired by a discusion with Boris Zbarsky
   and essentially just computes:

     h = 4.0/3.0 * tan ((angle_B - angle_A) / 4.0);

   without converting to polar coordinates.

   A different way to do this is covered in "Approximation of a cubic bezier
   curve by circular arcs and vice versa" by Alekas Riškus. However, the method
   presented there doesn't handle arcs with angles close to 0 because it
   divides by the perp dot product of the two angle vectors.

   NB: start and end must be normalized and cannot be elements of the output knots structure
*/
static double
_cairo_spline_arc (cairo_spline_knots_double_t *k, const cairo_vector_t *center, double radius, const cairo_vector_t *start, const cairo_vector_t *end)
{
    cairo_vector_t half = _cairo_vector_bisect (*start, *end);
    double coshalf = _cairo_vector_dot (*start, half);
    /* no need to use bisect since the angle between nstart and half is <= pi/2 < pi
     * and the vector doesn't need to be normalized */
    cairo_vector_t quart = _cairo_vector_add (*start, half);
    /* h = |start| |quart| sin(angle) / (|start| |quart| cos (angle) = tan (angle) */
    double h = (4./3.) * _cairo_vector_cross (*start, quart) / _cairo_vector_dot (*start, quart);

    k->a = _cairo_vector_scale (*start, radius);
    k->d = _cairo_vector_scale (*end, radius);
    k->b = _cairo_vector_add (k->a, _cairo_vector_scale(_cairo_vector_perpendicular (k->a), +h));
    k->c = _cairo_vector_add (k->d, _cairo_vector_scale(_cairo_vector_perpendicular (k->d), -h));
    _cairo_spline_translate (k, k, center);

    return radius * (7.-coshalf) * (2.+coshalf) * (2.+coshalf) / (27. * (1.+coshalf));
}

static void
_cairo_spline_normal_sweep (const cairo_spline_knots_double_t *k, double *sweep)
{
    cairo_vector_t hodograph[3];
    double dot;

    hodograph[0] = _cairo_vector_normalize (_cairo_vector_sub (k->b, k->a));
    hodograph[1] = _cairo_vector_normalize (_cairo_vector_sub (k->c, k->b));
    hodograph[2] = _cairo_vector_normalize (_cairo_vector_sub (k->d, k->c));

    /* the angle covered by each segment of the hodograph is < pi, thus its
     * tangent uniquely determines it */
    /* tan(a/2) = sin(a)/(1+cos(a)) */
    dot = _cairo_vector_dot (hodograph[0], hodograph[1]);
    if (dot != -1.0)
	sweep[0] = _cairo_vector_cross (hodograph[0], hodograph[1]) / (1.0 + dot);
    else
	sweep[0] = HUGE_VAL;

    dot = _cairo_vector_dot (hodograph[2], hodograph[1]);
    if (dot != -1.0)
	sweep[1] = _cairo_vector_cross (hodograph[2], hodograph[1]) / (1 + dot);
    else
	sweep[1] = HUGE_VAL;

}

/* split must be an array sufficient for the split curves (at most 3 splines, depending on
   the inflection points */
static int
_cairo_spline_split_inflection (cairo_spline_knots_double_t *split, cairo_bool_t *cusp, const cairo_spline_knots_double_t *k)
{
    double t_inflect[2];
    int r = 0, inflections = _cairo_spline_inflection_points (k, t_inflect);
    cairo_spline_knots_double_t remaining = *k;

    switch (inflections) {
    case 2:
	cusp[r] = _cairo_spline_de_casteljau (&remaining, split+r, &remaining, t_inflect[r]);
	r++;

	/* reparamaterize the second inflection point over the 'remaining' spline:
	 * (domain between inflection points)/(domain of remaining) */
	t_inflect[1] = (t_inflect[1] - t_inflect[0]) / (1.0 - t_inflect[0]);

	/* fall through */
    case 1:
	cusp[r] = _cairo_spline_de_casteljau (&remaining, split+r, &remaining, t_inflect[r]);
	r++;

	/* fall through */
    case 0:
	split[r++] = remaining;
	return r;

    default:
	ASSERT_NOT_REACHED;
	return 0;
    }
}

#endif /* CAIRO_SPLINE_PRIVATE_H */
