/* -*- Mode: c; tab-width: 8; c-basic-offset: 4; indent-tabs-mode: t; -*- */
#ifndef CAIRO_VECTOR_PRIVATE_H
#define CAIRO_VECTOR_PRIVATE_H

#include <math.h>

#if _XOPEN_SOURCE >= 600 || defined (_ISOC99_SOURCE)
#define ISFINITE(x) isfinite (x)
#else
#define ISFINITE(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

#include "cairoint.h"

typedef cairo_point_double_t cairo_vector_t;

static const cairo_vector_t _cairo_vector_zero = { 0, 0 };

static inline cairo_bool_t
_cairo_vector_equals (cairo_vector_t a, cairo_vector_t b)
{
    return a.x == b.x && a.y == b.y; /* TODO: fixed point precision? */
}

static inline cairo_bool_t
_cairo_vector_iszero (cairo_vector_t a)
{
    return _cairo_vector_equals (a, _cairo_vector_zero);
}

static inline cairo_bool_t
_cairo_vector_isfinite (cairo_vector_t a)
{
    return ISFINITE (a.x) && ISFINITE (a.y);
}

static inline double
_cairo_vector_cross (cairo_vector_t a, cairo_vector_t b)
{
    return a.x * b.y - a.y * b.x;
}

static inline double
_cairo_vector_dot (cairo_vector_t a, cairo_vector_t b)
{
    return a.x * b.x + a.y * b.y;
}

static inline double
_cairo_vector_length (cairo_vector_t v)
{
    return hypot (v.x, v.y);
}

static inline cairo_vector_t
_cairo_vector_opposite (cairo_vector_t p)
{
    cairo_vector_t r = { -p.x, -p.y };
    return r;
}

static inline cairo_vector_t
_cairo_vector_perpendicular (cairo_vector_t p)
{
    cairo_vector_t r = { -p.y, p.x};
    return r;
}

static inline cairo_vector_t
_cairo_vector_add (cairo_vector_t a, cairo_vector_t b)
{
    cairo_vector_t r = { a.x + b.x, a.y + b.y };
    return r;
}

static inline cairo_vector_t
_cairo_vector_sub (cairo_vector_t a, cairo_vector_t b)
{
    cairo_vector_t r = { a.x - b.x, a.y - b.y };
    return r;
}

static inline cairo_vector_t
_cairo_vector_scale (cairo_vector_t v, double scale)
{
    cairo_vector_t r = { v.x * scale, v.y * scale };
    return r;
}

static inline cairo_vector_t
_cairo_vector_normalize (cairo_vector_t v)
{
    assert (!_cairo_vector_iszero (v));
    return _cairo_vector_scale (v, 1 / _cairo_vector_length (v));
}

static inline cairo_vector_t
_cairo_vector_lerp (cairo_vector_t start, cairo_vector_t end, double t)
{
    return _cairo_vector_add (_cairo_vector_scale (start, 1-t), _cairo_vector_scale (end, t));
}

/* a and b must be normalized */
static cairo_vector_t
_cairo_vector_bisect_ccw (cairo_vector_t a, cairo_vector_t b)
{
    cairo_vector_t mid = _cairo_vector_add (a, b);

    if (_cairo_vector_iszero (mid))
	return _cairo_vector_perpendicular (a);

    mid = _cairo_vector_normalize (mid);
    /* If the angle between a and b is a reflex angle (> pi),
     * then the direction of the bisector is opposite to the
     * one we just computed. */
    if (_cairo_vector_cross (a, b) < 0)
	mid = _cairo_vector_opposite (mid);

    return mid;
}

/* a and b must be normalized */
static cairo_vector_t
_cairo_vector_bisect (cairo_vector_t a, cairo_vector_t b)
{
    cairo_vector_t mid = _cairo_vector_add (a, b);

    if (_cairo_vector_iszero (mid))
	return _cairo_vector_perpendicular (a);
    else 
	return _cairo_vector_normalize (mid);
}

#endif /* CAIRO_VECTOR_PRIVATE_H */
