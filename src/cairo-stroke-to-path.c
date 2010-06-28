/* -*- Mode: c; tab-width: 8; c-basic-offset: 4; indent-tabs-mode: t; -*- */
/* cairo - a vector graphics library with display and print output
 *
 * Copyright © 2009 Mozilla Foundation
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

/* Takes a path as input and produces a stroked path as output.
 * The resulting path may self intersect so it needs to be filled with
 * a non-zero winding rule. */

#define _BSD_SOURCE /* for hypot() */
#include "cairoint.h"
#include "cairo-path-fixed-private.h"
#include "cairo-spline-offset.h"

#define printf(a, ...)

typedef struct {
    cairo_path_fixed_t *path;
    cairo_matrix_t *ctm;

    /* we need to keep track of whether we've drawn anything yet
       so that we can do a move_to to the intial point of a curve
       before drawing it */
    cairo_bool_t has_current_point;
    cairo_int_status_t status;
} path_output_t;

static void line_to (path_output_t *path, double x, double y)
{
    if (path->ctm)
	cairo_matrix_transform_point (path->ctm, &x, &y);
    path->has_current_point = TRUE;

    if (!path->status)
	path->status = _cairo_path_fixed_line_to (path->path,
	                                          _cairo_fixed_from_double(x),
	                                          _cairo_fixed_from_double(y));
}

static void curve_to (path_output_t *path, double x0, double y0, double x1, double y1, double x2, double y2)
{
    if (path->ctm) {
	cairo_matrix_transform_point (path->ctm, &x0, &y0);
	cairo_matrix_transform_point (path->ctm, &x1, &y1);
	cairo_matrix_transform_point (path->ctm, &x2, &y2);
    }

    if (!path->status)
	path->status = _cairo_path_fixed_curve_to (path->path,
	    _cairo_fixed_from_double (x0), _cairo_fixed_from_double (y0),
	    _cairo_fixed_from_double (x1), _cairo_fixed_from_double (y1),
	    _cairo_fixed_from_double (x2), _cairo_fixed_from_double (y2));

}

static void close_path (path_output_t *path)
{
    if (!path->status)
	path->status = _cairo_path_fixed_close_path (path->path);
    if (!path->status)
	_cairo_path_fixed_new_sub_path (path->path);
}

static void
arc (path_output_t *path, double xc, double yc, double radius, vector_t a, vector_t b)
{
    /* find a vector that bisects the angle between a and b */
    vector_t mid_v = bisect (a, b);

    /* construct the arc using two curve segments */
    arc_segment (path, xc, yc, radius, a, mid_v);
    arc_segment (path, xc, yc, radius, mid_v, b);
}

static inline void
join_round (path_output_t *path, point_t center, vector_t a, vector_t b, double radius)
{
    /*
    int ccw = dot (perp (b), a) >= 0; // XXX: is this always true?
    yes, otherwise we have an interior angle.
    assert (ccw);
    */
    arc (path, center.x, center.y, radius, a, b);
}

static void
draw_circle (path_output_t *path, point_t center, double radius)
{
    vector_t a = {1, 0};
    vector_t b = {-1, 0};

    /* draw a line to the begining of the arc */
    line_to (path, center.x + radius * a.x, center.y + radius * a.y);

    arc (path, center.x, center.y, radius, a, b);
    arc (path, center.x, center.y, radius, b, a);
}

static point_t
c2p (cairo_matrix_t *ctm_inverse, cairo_point_t c)
{
    point_t p;
    p.x = _cairo_fixed_to_double(c.x);
    p.y = _cairo_fixed_to_double(c.y);
    if (ctm_inverse)
	cairo_matrix_transform_point (ctm_inverse, &p.x, &p.y);

    return p;
}

#if 1
static vector_t
compute_normal (point_t p0, point_t p1)
{
    vector_t n;
    /* u is parallel vector */
    double ux = p1.x - p0.x;
    double uy = p1.y - p0.y;

    /* we actually want the inverse square root here. then we won't have to divide */
    double ulen = hypot(ux, uy);
    /* cairo uses hypot which does sqrt(ux*ux + uy*uy) without the accuracy problems
       i.e. hypot is slower, but more accurate */
    assert (ulen != 0);
    // v is perpendicular *unit* vector
    n.x = -uy/ulen;
    n.y = ux/ulen;

    return n;
}
#else
#include <emmintrin.h>
/* a quick and dirty sse implementation. This is much lower
 * precision than the version above and not as fast as it
 * could be. */
static vector_t
compute_normal (point_t p0, point_t p1)
{
    vector_t n;
    /* u is parallel vector */
    float ux = p1.x - p0.x;
    float uy = p1.y - p0.y;


    __m128 mux = _mm_load_ss (&ux);
    __m128 muy = _mm_load_ss (&uy);

    __m128 mux2 = _mm_mul_ss(mux, mux);
    __m128 muy2 = _mm_mul_ss(muy, muy);

    __m128 rlen = _mm_rsqrt_ss(_mm_add_ss(muy2, mux2));
#if 1
    mux = _mm_mul_ss(mux, rlen);
    muy = _mm_sub_ps(_mm_setzero_ps(), _mm_mul_ss(muy, rlen)); //negate
    __m128 results = _mm_shuffle_ps (mux, muy, _MM_SHUFFLE(0, 0, 0, 0)); // get bottom words
    results = _mm_shuffle_ps (results, results, _MM_SHUFFLE(0, 0, 0, 2)); // x and y
#else
    rlen = _mm_shuffle_ps (rlen, rlen, _MM_SHUFFLE(0, 0, 0, 0)); // get bottom words
    muy = _mm_sub_ps(_mm_setzero_ps(), muy);
    __m128 results = _mm_shuffle_ps (mux, muy, _MM_SHUFFLE(0, 0, 0, 0)); // get bottom words
    results = _mm_shuffle_ps (results, results, _MM_SHUFFLE(0, 0, 0, 2)); // x and y
    results = _mm_mul_ps(mux, rlen);
#endif


    __m128d result = _mm_cvtps_pd (results);
    _mm_storeu_pd(&n.x, result);
#if 0
    __m128 muxy = _mm_shuffle_ps (mux, muy, _MM_SHUFFLE(0, 0, 0, 0)); // get bottom words
    /* we actually want the inverse square root here. then we won't have to divide */
    muxy = _mm_mul_ps(muxy, muxy);
    muxy = _mm_shuffle_ps (mux, muy, _MM_SHUFFLE(0, 0, 0, 0)); // get bottom words
#endif
    return n;
}
#endif



static void offset_curve_to(void *closure, knots_t c)
{
    path_output_t *path = closure;
    ensure_finite(c);
    //XXX: it would be nice if we didn't have to check this for every point
    if (!path->has_current_point)
	move_to(path, c.a.x, c.a.y);
    curve_to(path, c.b.x, c.b.y, c.c.x, c.c.y, c.d.x, c.d.y);
}

static void
draw_offset_curve (path_output_t *path,
	point_t a, point_t b, point_t c, point_t d, double offset)
{
    knots_t curve = {a, b, c, d};
    /* we want greater tolerance when the lines we are stroking are narrower
     * because any deviations will be more noticeable.
     * Fortunately, the narrow the lines are the easier it is to approximate
     * the offset curve. It could be argued that this is a bad choice because
     * one could draw two large but slightly different stroke widths on top of
     * each other and the result wouldn't necessarily match. However,
     * I'm not sure if tha's actually a problem. */
    double tolerance = offset * 2 / 500.; /* the choice of 500 is pretty arbitrary */

    curve_offset(curve, offset, tolerance, offset_curve_to, path);
}

/* Finds the intersection of two lines each defined by a point and a normal.
   From "Example 2: Find the interesection of two lines" of
   "The Pleasures of "Perp Dot" Products"
   F. S. Hill, Jr. */
static point_t
line_intersection (point_t A, vector_t a_perp, point_t B, vector_t b_perp)
{
    point_t intersection;
    double denom;
    double t;

    vector_t a = unperp(a_perp);
    vector_t c = {B.x - A.x, B.y - A.y};
    denom = dot(b_perp, a);
    /* CHECK ME!!! */
    if (denom == 0.0) {
	assert(0 && "TROUBLE");
    }

    t = dot (b_perp, c)/denom;

    intersection.x = A.x+t*(a.x);
    intersection.y = A.y+t*(a.y);

    return intersection;
}

static cairo_bool_t
is_interior_angle (vector_t a, vector_t b) {
    /* angles of 180 and 0 degress will evaluate to 0, however
     * we to treat 180 as an interior angle and 180 as an exterior angle */
    return dot (perp (a), b) > 0 ||
	(a.x == b.x && a.y == b.y); /* 0 degrees is interior */
}

static void
join_segment_line (path_output_t *path, cairo_stroke_style_t *style, point_t p, vector_t s1_normal, vector_t s2_normal)
{
    double miter_limit = style->miter_limit;
    cairo_line_join_t join = style->line_join;
    double offset = style->line_width/2;

    point_t start = {p.x + s1_normal.x*offset, p.y + s1_normal.y*offset};
    point_t end   = {p.x + s2_normal.x*offset, p.y + s2_normal.y*offset};

    if (is_interior_angle (s1_normal, s2_normal)) {
	line_to(path, start.x, start.y);
	//XXX: while you wouldn't think this point is necessary, it is when we the segments we are joining
	//are shorter than the line width
	/* qt does a check to avoid adding this additional point when possible */
	line_to(path, p.x, p.y);

	line_to(path, end.x, end.y);
    } else {
	switch (join) {
	    case CAIRO_LINE_JOIN_ROUND:
		line_to(path, start.x, start.y);
		join_round(path, p, s1_normal, s2_normal, offset);
		break;
	    case CAIRO_LINE_JOIN_MITER:
		{
		    double in_dot_out = -s1_normal.x*s2_normal.x + -s1_normal.y*s2_normal.y;
		    // XXX: we are going to have an extra colinear segment here */
		    if (2 <= miter_limit*miter_limit * (1 - in_dot_out)) {
			point_t intersection = line_intersection(start, s1_normal, end, s2_normal);
			line_to(path, intersection.x, intersection.y);
			break;
		    }
		}
		// fall through
	    case CAIRO_LINE_JOIN_BEVEL:
		line_to(path, start.x, start.y);
		line_to(path, end.x, end.y);
		break;
	}
    }
}


static void
join_segment_curve (path_output_t *path, cairo_stroke_style_t *style, point_t p, vector_t s1_normal, vector_t s2_normal)
{
    cairo_line_join_t join = style->line_join;
    double offset = style->line_width/2;
    double miter_limit = style->miter_limit;

    point_t start = {p.x + s1_normal.x*offset, p.y + s1_normal.y*offset};
    point_t end   = {p.x + s2_normal.x*offset, p.y + s2_normal.y*offset};
    if (is_interior_angle (s1_normal, s2_normal)) {
	//XXX: while you wouldn't think this point is necessary, it is when we the segments we are joining
	//are shorter than the line width
	/* qt does a check to avoid adding this additional point when possible */
	line_to(path, p.x, p.y);
	line_to(path, end.x, end.y);
    } else {
	switch (join) {
	    case CAIRO_LINE_JOIN_ROUND:
		join_round(path, p, s1_normal, s2_normal, offset);
		break;
	    case CAIRO_LINE_JOIN_MITER:
	    {
		double in_dot_out = -s1_normal.x*s2_normal.x + -s1_normal.y*s2_normal.y;
		// XXX: we are going to have an extra colinear segment here */
		if (2 <= miter_limit*miter_limit * (1 - in_dot_out)) {
		    point_t intersection = line_intersection(start, s1_normal, end, s2_normal);
		    line_to(path, intersection.x, intersection.y);
		    break;
		}
	    }
		// fall through
	    case CAIRO_LINE_JOIN_BEVEL:
		line_to(path, end.x, end.y);
		break;
	}
    }
}

/* we have some different options as to what the semantics of cap_line should be:
   - here's how skia's works:
     - it uses a virtual function call instead of a switch statement.
     - it assumes that the caller has added the start point of the cap
     - here's an interesting bit: the square capper has different behaviour
       depending on whether the last segment is a curve or not.
       If it is a line it will adjust the last point of the path
       and only set the outside point of the cap. This avoids the
       the extra points when capping a square line.
  - qt and mesa's openvg state tracker do this cute thing
    where they treat caps as joins. i.e a butt cap maps
    to a miter join, a round cap to a round join and
    a square cap is special cased.
    3 line segments are always emitted for a square cap
    it seems the code assumes that the first point is implied?

    I don't think the similarities of caps and joins outweigh the differences. Further
    I think keeping them separate makes it easier to understand.
*/
static void
cap_line (path_output_t *path, cairo_stroke_style_t *style, point_t end, vector_t normal)
{
    cairo_line_cap_t cap = style->line_cap;
    double offset = style->line_width/2;

    /* This function will add additional unnecessary lines when it is called after a curve
     * segment. I'm not sure it's worth trying to eliminate them */
    switch (cap) {
	case CAIRO_LINE_CAP_ROUND:
	    {
		line_to (path, end.x + offset*normal.x, end.y + offset*normal.y);
		arc(path, end.x, end.y, offset, normal, flip(normal));
		break;
	    }
	case CAIRO_LINE_CAP_SQUARE:
	    {
		// parallel vector
		double vx = normal.y;
		double vy = -normal.x;
		// move the end point down the line by offset
		end.x += vx*offset;
		end.y += vy*offset;

		// offset the end point of last_line to draw the cap at the end
		line_to(path, end.x + offset*normal.x, end.y + offset*normal.y);
		line_to(path, end.x - offset*normal.x, end.y - offset*normal.y);
		break;
	    }
	case CAIRO_LINE_CAP_BUTT:
	    {
		line_to (path, end.x + offset*normal.x, end.y + offset*normal.y);
		// offset the end point of last_line to draw the cap at the end
		line_to(path, end.x - offset*normal.x, end.y - offset*normal.y);
		break;
	    }
    }
}

#define cairo_path_head(path__) (&(path__)->buf.base)
#define cairo_path_tail(path__) cairo_path_buf_prev (cairo_path_head (path__))

#define cairo_path_buf_next(pos__) \
    cairo_list_entry ((pos__)->link.next, cairo_path_buf_t, link)
#define cairo_path_buf_prev(pos__) \
    cairo_list_entry ((pos__)->link.prev, cairo_path_buf_t, link)

// index is the first segment we need to deal with
typedef struct {
    const cairo_path_buf_t *buf;
    cairo_point_t *points;
    unsigned int op_index;
    const cairo_path_buf_t *first;
} path_ittr;

static const uint8_t num_args[] = {
	1, /* cairo_path_move_to */
	1, /* cairo_path_op_line_to */
	3, /* cairo_path_op_curve_to */
	0, /* cairo_path_op_close_path */
    };

static inline path_ittr
ittr_next (path_ittr index)
{
    index.points += num_args[(int)index.buf->op[index.op_index]];
    index.op_index++;
    if (index.op_index == index.buf->num_ops) {
	index.buf = cairo_path_buf_next (index.buf);
#ifdef ITTR_DEBUG
	if (index.buf == index.first) {
	    assert(0);
	    index.buf = NULL;
	    index.points = NULL;
	    index.op_index = 0xbaadf00d;
	}
#endif
	index.op_index = 0;
	index.points = index.buf->points;
    }
    return index;
}

static inline path_ittr
ittr_prev (path_ittr index)
{
    if (index.op_index == 0) {
	index.buf = cairo_path_buf_prev (index.buf);
	index.op_index = index.buf->num_ops-1;
	index.points = index.buf->points + (index.buf->num_points);
    } else {
	index.op_index--;
    }
    index.points -= num_args[(int)index.buf->op[index.op_index]];
    return index;
}

static inline cairo_path_op_t
ittr_op (path_ittr i)
{
    return i.buf->op[i.op_index];
}

static cairo_bool_t
ittr_done(path_ittr i)
{
    return i.buf == NULL;
}

static cairo_bool_t
ittr_last (path_ittr i)
{
    return i.op_index+1 == i.buf->num_ops && cairo_path_buf_next(i.buf) == i.first;
}

static cairo_bool_t
ittr_eq (path_ittr a, path_ittr b)
{
    return a.buf == b.buf && a.op_index == b.op_index;
}


static void
draw_degenerate_path (path_output_t *path_out, cairo_stroke_style_t *style, point_t point)
{
    if (style->line_cap == CAIRO_LINE_CAP_ROUND) {
	draw_circle(path_out, point, style->line_width/2);
    }
}

/* it might be interesting to see if we can just skip degenerate segments... */

static cairo_status_t
_cairo_subpath_stroke_to_path (path_ittr *ittr,
			    cairo_stroke_style_t *stroke_style,
			    path_output_t *path_out,
			    cairo_matrix_t *ctm_inverse)
{
    double offset = stroke_style->line_width/2;
    cairo_status_t status = CAIRO_STATUS_SUCCESS;
    path_ittr i = *ittr;
    cairo_matrix_t *ctm_i = ctm_inverse;

    point_t begin = c2p (ctm_i, i.points[0]);
    point_t end_point = begin;

    cairo_path_op_t op;
    cairo_point_t *points;
    vector_t end_normal;
    vector_t initial_normal;
    point_t end;

    int op_count = 0;
    cairo_bool_t closed_path = FALSE;
    point_t last_point;
    vector_t closed_normal;

    assert(ittr_op(i) == CAIRO_PATH_OP_MOVE_TO);

    /* handle beginging of the line: we want to get an intial normal to work with
     * before we start drawing */
    while (1) {
	i = ittr_next (i);
	op = ittr_op (i);
	points = i.points;

	if (op == CAIRO_PATH_OP_MOVE_TO) {
	    *ittr = i;
	    draw_degenerate_path(path_out, stroke_style, begin);
	    return CAIRO_STATUS_SUCCESS;
	}

	if (op == CAIRO_PATH_OP_CLOSE_PATH) {
	    *ittr = ittr_next (i);
	    draw_degenerate_path(path_out, stroke_style, begin);
	    return CAIRO_STATUS_SUCCESS;
	}

	if (op == CAIRO_PATH_OP_LINE_TO) {
	    point_t p = c2p(ctm_i, points[0]);
	    //XXX: it would be better if we compared points in fixed point...
	    if (!point_eq(begin, p)) {
		initial_normal = compute_normal(begin, p);
		break;
	    }
	} else {
	    point_t p0 = c2p(ctm_i, points[0]);
	    point_t p1 = c2p(ctm_i, points[1]);
	    point_t p2 = c2p(ctm_i, points[2]);
	    assert(ittr_op (i) == CAIRO_PATH_OP_CURVE_TO);
	    //XXX: maybe use the helper function?
	    if (!point_eq (begin, p0)) {
		initial_normal = compute_normal (begin, p0);
		break;
	    } else {
		if (!point_eq (begin, p1)) {
			initial_normal = compute_normal (begin, p1);
			break;
		} else {
		    if (!point_eq (begin, p2)) {
			initial_normal = compute_normal (begin, p2);
			break;
		    }
		}
	    }
	}

	if (ittr_last (i)) {
	    /* if we had a normal for the last op we would have broken out already */
	    printf("degen\n");
	    *ittr = i;
	    draw_degenerate_path (path_out, stroke_style, begin);
	    return CAIRO_STATUS_SUCCESS;
	}
    }

    end_normal = initial_normal;

    /* end_normal needs to be the normal between i.points[] and the preceeding point */
    while (!ittr_last (i)) {
	    op_count++;

	    i = ittr_next (i);
	    cairo_path_op_t next_op = ittr_op (i);
	    cairo_point_t *next_points = i.points;

	    switch (next_op) {
	    case CAIRO_PATH_OP_MOVE_TO:
		*ittr = i;
		i = ittr_prev (i);
		op_count--;
		// we are done with the current subpath
		// cap it and return to the begining
		//assert(0);
		goto done;
	    case CAIRO_PATH_OP_LINE_TO:
	    case CAIRO_PATH_OP_CURVE_TO:
		// get first point of next_op
		end = c2p (ctm_i, next_points[0]);
		break;
	    case CAIRO_PATH_OP_CLOSE_PATH:
		// close the path
		// XXX: relies on the fact that CLOSE_PATH is never the last path_op (it is always followed by an OP_MOVE)
		closed_path = TRUE;
		// the first point of the first op
		end = c2p (ctm_i, ittr->points[0]);
		break;
	    default:
		ASSERT_NOT_REACHED;
	    }
	    if (unlikely (status))
		return status;

	    switch (op) {
	    case CAIRO_PATH_OP_CURVE_TO:
		{
		    point_t p0 = c2p (ctm_i, points[0]);
		    point_t p1 = c2p (ctm_i, points[1]);
		    point_t p2 = c2p (ctm_i, points[2]);

		    vector_t n2 = end_normal;
		    if (!point_eq (p0, p1))
			n2 = compute_normal (p0, p1);
		    if (!point_eq (p1, p2))
			n2 = compute_normal (p1, p2);
		    end_normal = n2;

		    draw_offset_curve (path_out, begin,
			    p0,
			    p1,
			    p2,
			    offset);

		    if (!point_eq (p2, end)) {
			vector_t n3 = compute_normal (p2, end);
			join_segment_curve (path_out, stroke_style, p2, n2, n3);
			end_normal = n3;
		    }

		    begin = p2;
		    end_point = begin;
		}
		break;
	    case CAIRO_PATH_OP_LINE_TO:
		{
		    point_t mid = c2p (ctm_i, points[0]);
		    if (!point_eq (mid, end)) {
			vector_t n1 = end_normal;
			vector_t n2 = compute_normal (mid, end);
			join_segment_line (path_out, stroke_style, mid, n1, n2);
			end_normal = n2;
			begin = mid;
			end_point = mid;
		    }
		}
		break;
	    default:
		ASSERT_NOT_REACHED;
	    }
	    if (closed_path)
		break;
	    op = next_op;
	    points = next_points;
    };

    if (!closed_path)
	*ittr = i;
done:
    // op should be set the last op in the path
    if (closed_path) {
	point_t first_point = c2p (ctm_i, ittr->points[0]);

	assert(ittr_op (*ittr) == CAIRO_PATH_OP_MOVE_TO);
	closed_normal = end_normal;
	last_point = end_point;
	join_segment_line (path_out, stroke_style, first_point, end_normal, initial_normal);

	// XXX: is this right?
	begin = first_point;
	*ittr = ittr_next (i);
	path_out->has_current_point = FALSE;
	close_path (path_out);
    } else {
	if (op == CAIRO_PATH_OP_CURVE_TO) {
	    point_t p0 = c2p (ctm_i, points[0]);
	    point_t p1 = c2p (ctm_i, points[1]);
	    point_t p2 = c2p (ctm_i, points[2]);

	    /* draw capped curve segment */
	    draw_offset_curve (path_out, end_point, p0, p1, p2, offset);
	    cap_line (path_out, stroke_style, p2, curve_normal (end_normal, end_point, p0, p1, p2));
	    draw_offset_curve (path_out, p2, p1, p0, end_point, offset);

	    begin = p0;
	} else if (op == CAIRO_PATH_OP_LINE_TO) {
	    begin = c2p (ctm_i, points[0]);
	    cap_line (path_out, stroke_style, begin, end_normal);
	} else {
	    assert(0);
	}
    }

    end_normal = flip (end_normal);

    // pre-conditions
    i = ittr_prev(i);
    op = ittr_op(i);
    points = i.points;

    /* traverse backwards to the begining drawing
     * the other side of the stroked path */
    while (op_count) {
	i = ittr_prev (i);
	cairo_path_op_t next_op = ittr_op (i);
	cairo_point_t *next_points = i.points;

	switch (next_op) {
	    case CAIRO_PATH_OP_MOVE_TO:
		end = c2p (ctm_i, next_points[0]);
		break;
	    case CAIRO_PATH_OP_LINE_TO:
		end = c2p (ctm_i, next_points[0]);
		break;
	    case CAIRO_PATH_OP_CURVE_TO:
		// get first point of next_op
		end = c2p (ctm_i, next_points[2]);
		break;
	    case CAIRO_PATH_OP_CLOSE_PATH:
		// close the path
		assert(0);
		break;
	    default:
		ASSERT_NOT_REACHED;
	}
	if (unlikely (status))
	    return status;

	switch (op) {
	    case CAIRO_PATH_OP_CURVE_TO:
		{
		    point_t p2 = c2p (ctm_i, points[2]);
		    point_t p1 = c2p (ctm_i, points[1]);
		    point_t p0 = c2p (ctm_i, points[0]);
		    if (!point_eq (p2, p1)) {
			vector_t n0 = compute_normal (p2, p1);
			// XXX whether we use join_segment_line or join_segment_curve here depends
			// on whether we are coming from a line or from a curve. Does it matter?
			// this used to be join_segment_curve but it was wrong for the following:
			//         cairo_move_to(cr, 550, 400);
			//	   cairo_curve_to(cr, 531.687500, 469.878906, 415.203125, 506.437500, 247.179688, 524.625000);
			//         cairo_line_to(cr, 257.945312, 624.046875);
			join_segment_line (path_out, stroke_style, p2, end_normal, n0);
			end_normal = n0;
		    }
		    end_normal = curve_normal (end_normal, p2, p1, p0, end);
		    draw_offset_curve (path_out,
			    p2,
			    p1,
			    p0,
			    end,
			    offset);
		    //XXX: this feels suspect...
		    begin = p0;
		    end_point = end;
		}
		break;
	    case CAIRO_PATH_OP_LINE_TO:
		{
		    point_t mid = c2p (ctm_i, points[0]);
		    vector_t n1 = end_normal;
		    if (!point_eq (mid, end)) {
			vector_t n2 = compute_normal (mid, end);
			join_segment_line (path_out, stroke_style, mid, n1, n2);
			end_normal = n2;
			begin = mid;
			end_point = end;
		    }
		}
		break;
	    default:
		ASSERT_NOT_REACHED;
	}

	points = next_points;
	op = next_op;
	op_count--;

    }

    if (closed_path) {
	vector_t n2 = end_normal;
	printf("finish closed\n");
	assert(ittr_op (i) == CAIRO_PATH_OP_MOVE_TO);
	if (!point_eq (end_point, last_point)) {
	    n2 = compute_normal (end_point, last_point);
	    join_segment_line (path_out, stroke_style, c2p (ctm_i, i.points[0]), end_normal, n2);
	}

	join_segment_line (path_out, stroke_style, last_point, n2, flip (closed_normal));
    } else {
	cap_line (path_out, stroke_style, end_point, end_normal);
    }
    close_path (path_out);

    return CAIRO_STATUS_SUCCESS;
}

typedef struct {
        point_t point;
        vector_t normal;
} dash_point_t;

static double
distance (point_t p0, point_t p1)
{
    double ux = p1.x - p0.x;
    double uy = p1.y - p0.y;
    return hypot(ux, uy);
}

/* i points to the last op we've just finished */
static dash_point_t
skip_subpath (dash_point_t start_point, path_output_t *path_out, path_ittr *i, path_ittr subpath_start, path_ittr *subpath_end, cairo_matrix_t *ctm_i, cairo_bool_t *closed, double length)
{
    //XXX: review the naming of points to make sure it makes sense and is consistent
    path_ittr index = *i;
    path_ittr next_index;
    dash_point_t end_point = start_point;
    point_t next_point;
    cairo_bool_t done = FALSE;

    assert(!ittr_last(index));

    next_index = ittr_next(index);

    printf("skip_subpath: %d\n", index.op_index);

    if (ittr_op(next_index) == CAIRO_PATH_OP_CLOSE_PATH) {
	*closed = TRUE;
	next_index = subpath_start;
    }

    assert(ittr_op(next_index) == CAIRO_PATH_OP_LINE_TO || ittr_eq(next_index, subpath_start));
    next_point = c2p(ctm_i, next_index.points[0]);

    while (length > 0) { //XXX: is this condition even needed?
	double d = distance(end_point.point, next_point);
	assert(ittr_op(next_index) == CAIRO_PATH_OP_LINE_TO || ittr_op(next_index) == CAIRO_PATH_OP_MOVE_TO);
	printf("d(%d): %f of %f\n", index.op_index, d, length);
	if (d >= length) {
	    /* move the end_point a distance 'length' along the line with the normal end_point.normal */
	    end_point.point.x += -end_point.normal.y*-length;
	    end_point.point.y += end_point.normal.x*-length;
	    length = 0;
	    printf("skippath break\n");
	    break;
	} else {
	    end_point.point = next_point;

	    if (ittr_eq(next_index, subpath_start)) {
	        /* we've passed the begining; we are done */
		done = TRUE;
		break;
	    }

	    if (ittr_last(next_index)) {
		/* we're at the end of the path; we are done */
		*subpath_end = next_index;
		done = TRUE;
		break;
	    }

	    index = next_index;
	    next_index = ittr_next(next_index);

	    if (ittr_op(next_index) == CAIRO_PATH_OP_CLOSE_PATH) {
		/* we have a closed path so wrap around toward the begining (subpath_start) */
		*closed = TRUE;
		*subpath_end = ittr_next(next_index);
		next_index = subpath_start;
	    } else if (ittr_op(next_index) == CAIRO_PATH_OP_MOVE_TO) {
		*subpath_end = next_index;
		done = TRUE;
		break;
	    }

	    next_point  = c2p(ctm_i, next_index.points[0]);
	    if (!point_eq(end_point.point, next_point)) {
		/* XXX: we might be able to avoid computing the normal for every point we vist
		 * and instead only compute it for the last segment */
		end_point.normal = compute_normal(end_point.point, next_point);
	    }

	    printf("next %f %f\n", d, length);
	    length -= d;
	}
    }

    *i = index;

    if (done)
	i->buf = NULL;

    return end_point;
}

static dash_point_t
dash_subpath (dash_point_t start_point, path_output_t *path_out, path_ittr *index, path_ittr subpath_start, path_ittr *subpath_done, cairo_matrix_t *ctm_i, cairo_stroke_style_t *style, double length, int begin_on, cairo_bool_t *closed, double start_length)
{
    dash_point_t end_point;
    point_t first = start_point.point;
    path_ittr i = *index;
    path_ittr subpath_end = {};
    int closed_path = 0;
    int done = 0;
    vector_t n1 = start_point.normal;
    path_ittr start_index = i;
    path_ittr next_i;
    point_t p1;
    // investigate more normal reuse.
    // ideally we would only compute the normal for each segment once
    // deal with short paths (count < 3)
    assert(!ittr_last(i));
    printf("dash_subpath(%f): %d: %f %f\n", length, i.op_index, first.x, first.y);

    next_i = ittr_next(i);

    if (ittr_op(next_i) == CAIRO_PATH_OP_CLOSE_PATH) {
	subpath_end = i;
	*subpath_done = ittr_next(next_i);
	next_i = subpath_start;
	*closed = TRUE;
    }

    assert(ittr_op(next_i) == CAIRO_PATH_OP_LINE_TO || ittr_eq(next_i, subpath_start));
    //printf("off: %f %f, %f %f\n", end_point.point.x, end_point.point.y, path[index].x, path[index].y);
    p1 = c2p(ctm_i, next_i.points[0]);
    do {
	double d;
	vector_t n2;
	point_t p2;
	assert(ittr_op(next_i) == CAIRO_PATH_OP_LINE_TO || ittr_eq(next_i, subpath_start));
	d = distance (first, p1);
	if (d >= length) {
	    end_point.point.x = first.x + length*n1.y;
	    end_point.point.y = first.y + length*-n1.x;
	    end_point.normal = n1;
	    break;
	}

	if (ittr_last (next_i)) {
	    /* we're at the end so we must be done */
	    done = 1;
	    *subpath_done = next_i;
	    end_point.point = c2p(ctm_i, next_i.points[0]);
	    end_point.normal = n1;
	    break;
	}

	/* We don't have enough length in the current segment so get the next one
	 * This means advancing i and next_i */

	if (ittr_eq(next_i, subpath_start)) {
	    /* we're already at the beginning of the path so we don't have anything left to dash */
	    printf("done\n");
	    //XXX done the path
	    done = 1;
	    if (!begin_on) {
		/* we didn't skip the initial segment; so we're all done */
		length = 0;
		end_point.point = c2p(ctm_i, next_i.points[0]);
		end_point.normal = n1;
		break;
	    } else {
		/* the initial segment started on, so join with it and draw it */
		length = start_length;
		i = next_i;
		next_i = ittr_next(next_i);
	    }
	} else {
	    i = next_i;
	    next_i = ittr_next(next_i);

	    if (ittr_op(next_i) == CAIRO_PATH_OP_CLOSE_PATH) {
		subpath_end = i;
		*subpath_done = ittr_next(next_i);
		next_i = subpath_start;
		*closed = TRUE;
	    } else if (ittr_op(next_i) == CAIRO_PATH_OP_MOVE_TO) {
		done = 1;
		end_point.point = c2p(ctm_i, i.points[0]);
		end_point.normal = n1;
		*subpath_done = next_i;
		//XXX: it would be good if we could avoid doing this ittr_prev()
		i = ittr_prev(i);
		printf("dash: segment end break\n");
		break;
	    }
	    length -= d;
	}

	p2 = c2p(ctm_i, next_i.points[0]);
	if (!point_eq(p1, p2)) {
	    n2 = compute_normal(p1, p2);
	    end_point.normal = n2;
	    end_point.point = p2;
	}

	if (ittr_eq(next_i, start_index)) {
	    // we've reached the place where we started so we need to loop the path instead of
	    // capping it
	    done = 1;
	    closed_path = 1;
	    /*
	    gcc complains that end_point can be used uninitialzed in this function because we don't set
	    it here when we break out of the loop. I don't think we can get to here without setting end_point
	    and here's why:

	    p2 == start_index point
	    p1 must have the same value as p2

	    the only time that next_i can equal start_index is when start_index == subpath_start
	    this means that the whole path needs to be here. However, unless the path is degenerate
	    end_point will be set. But no degenerate paths will pass into here because they will be
	    found else where.

	    Can we make this more explict?
	    */
	    printf("closed path\n");
	    // XXX: do we need to set end_point here? normally we would get an end_point from the code above.
	    // but if we have a degenerate path we won't
	    break;
	}

	// XXX: having to check this twice sort of sucks
	if (!point_eq(p1, p2)) {
	    join_segment_line(path_out, style, p1, n1, n2);
	    n1 = n2;
	    first = p1;
	    p1 = p2;
	}
    } while (1);

    /* at this point end_point must be initialized and 'i' must point to the point before the end_point */

    *index = i;

    point_t closed_p1;
    point_t closed_p2;
    vector_t closed_n1;
    vector_t closed_n2;
    vector_t closed_n3;
    if (closed_path) {
	    // join the line back with itself
	    closed_p1 = c2p (ctm_i, i.points[0]);
	    closed_p2 = end_point.point;
	    closed_n1 = n1;
	    /* CHECK ME!!! */
	    closed_n2 = compute_normal (closed_p1, closed_p2);
	    closed_n3 = compute_normal (closed_p2, c2p(ctm_i, ittr_next (next_i).points[0]));
	    join_segment_line (path_out, style, closed_p1, closed_n1, closed_n2);
	    join_segment_line (path_out, style, closed_p2, closed_n2, closed_n3);
	    close_path (path_out);
    } else {
	// draw cap
	cap_line (path_out, style, end_point.point, end_point.normal);
    }

    start_point.normal = flip (start_point.normal);

    /* draw the other side of the path, returning to the begining */
    n1 = flip (end_point.normal);
    p1 = c2p (ctm_i, i.points[0]);
    path_ittr prev_i = i;
    while (!ittr_eq (prev_i, start_index)) {
	i = prev_i;
	if (ittr_eq (prev_i, subpath_start)) {
	    prev_i = subpath_end;
	} else {
	    prev_i = ittr_prev (prev_i);
	}
	point_t p2 = c2p (ctm_i, prev_i.points[0]);
	if (!point_eq (p1, p2)) {
	    vector_t n2 = compute_normal (p1, p2);
	    join_segment_line (path_out, style, p1, n1, n2);

	    n1 = n2;
	    p1 = p2;
	}
    }

    /* finish up */

    if (closed_path) {
	join_segment_line (path_out, style, closed_p2, flip (closed_n3), flip (closed_n2));
	join_segment_line (path_out, style, closed_p1, flip (closed_n2), flip (closed_n1));
	close_path (path_out);
    } else {
	/* draw cap */
	cap_line(path_out, style, start_point.point, start_point.normal);
	close_path (path_out);
    }

    if (done)
	index->buf = NULL;

    return end_point;

}

static cairo_status_t
_cairo_subpath_dash_to_path (path_ittr *ittr,
			    cairo_stroke_style_t *style,
			    path_output_t *path_out,
			    cairo_matrix_t *ctm_i)
{

    cairo_bool_t on = TRUE, begin_on, closed = FALSE;

    double length, start_length;

    double *dash = style->dash;
    double dash_init = style->dash_offset;
    int dash_len = style->num_dashes;
    int dash_state = 0;
    path_ittr index = *ittr;
    path_ittr subpath_end = {};

    dash_point_t start_point;
    dash_point_t end_point;
    path_ittr subpath_start = index;

    assert(ittr_op(index) == CAIRO_PATH_OP_MOVE_TO);

    ///XXX: this code is different from cairo, figure out which is better
    /* intiailize dash state */
    length = dash[dash_state++ % dash_len];

    /* decrease dash_init until we are in the correct dash_state */
    /* the intervals are closed */
    /* we make sure that we don't skip over 0 length segments */
    while (length > 0 && dash_init - length >= 0) {
	dash_init -= length;
	on = !on;
	length = dash[dash_state++ % dash_len];
    }

    begin_on = on;
    length -= dash_init;
    printf("subpath starts %s with length of %f\n", on ? "on" : "off", length);

    start_length = length;

    start_point.point = c2p(ctm_i, index.points[0]);
    //XXX: to share this with the non dashed code we need to deal with dash on/off
    while (1) {
	if (ittr_op(ittr_next(index)) == CAIRO_PATH_OP_LINE_TO) {
	    point_t p = c2p(ctm_i, ittr_next(index).points[0]);
	    if (!point_eq(start_point.point, p)) {
		/* we have a normal so we can start the regular stroking process */
		start_point.normal = compute_normal(start_point.point, c2p(ctm_i, ittr_next(index).points[0]));
		break;
	    }
	}
	index = ittr_next(index);

	if (ittr_op(index) == CAIRO_PATH_OP_MOVE_TO) {
	    *ittr = index;
	    if (on)
		draw_degenerate_path (path_out, style, start_point.point);
	    return CAIRO_STATUS_SUCCESS;
	}

	if (ittr_op(index) == CAIRO_PATH_OP_CLOSE_PATH) {
	    *ittr = ittr_next(index);
	    if (on)
		draw_degenerate_path (path_out, style, start_point.point);
	    return CAIRO_STATUS_SUCCESS;
	}

	if (ittr_last(index)) {
	    *ittr = index;
	    if (on)
		draw_degenerate_path (path_out, style, start_point.point);
	    return CAIRO_STATUS_SUCCESS;
	}
    }

    end_point = start_point;

    /* the first segment is special because we might need to join with it. However, we won't know
     * what to do until we reach the end of the subpath
     * Different renderers behave differently here:
     *   No Join: skia, qt, mesa openvg, coregraphics
     *   Join: opera, cairo, ghostscript, adobe */
    if (begin_on) {
	on = !on;
	/* if the first dash segment is larger than the entire line we'll skip the whole thing */
	end_point = skip_subpath (end_point, path_out, &index, subpath_start, &subpath_end, ctm_i, &closed, length);
	length = dash[dash_state++ % dash_len];
	printf("done skip: %d\n", index.op_index);
	if (ittr_done(index)) {
	    printf("reset\n");
	    // flip to back on so that we don't skip everything again
	    on = !on;
	    index = *ittr;
	    length = start_length;
	    end_point = start_point;
	}
    }

    do {
	printf("index: %d\n", index.op_index);
	if (on) {
	    end_point = dash_subpath (end_point, path_out, &index, subpath_start, &subpath_end, ctm_i, style, length, begin_on, &closed, start_length);
	} else {
	    end_point = skip_subpath (end_point, path_out, &index, subpath_start, &subpath_end, ctm_i, &closed, length);
	}
	on = !on;
	length = dash[dash_state++ % dash_len];

	// cairo does an implicit move_to after a close_path so we need to start a new sub path
	_cairo_path_fixed_new_sub_path (path_out->path);
    } while (!ittr_done(index));

    *ittr = index;
    // draw the initial segment if we haven't yet.
    if ((!closed && begin_on) || (closed && on)) {
	printf("draw initial %f\n", start_length);
	index = subpath_start;
	end_point = start_point;
	end_point = dash_subpath (end_point, path_out, &index, subpath_start, &subpath_start, ctm_i, style, start_length, begin_on, &closed, start_length);
    }
    *ittr = subpath_end;
    return CAIRO_STATUS_SUCCESS;
}

static cairo_status_t
_cairo_path_dash_to_path (const cairo_path_fixed_t		*path,
			    cairo_stroke_style_t *stroke_style,
			     cairo_path_fixed_t         *path_out,
			     cairo_matrix_t *ctm,
			     cairo_matrix_t *ctm_inverse)
{
    path_ittr i;
    path_output_t output_path;

    output_path.status = CAIRO_STATUS_SUCCESS;
    output_path.ctm = ctm;
    output_path.path = path_out;

    i.buf = i.first = cairo_path_head(path);
    i.op_index = 0;
    i.points = i.buf->points;

    while (!ittr_last(i)) {
	_cairo_subpath_dash_to_path(&i, stroke_style, &output_path, ctm_inverse);
    }

    return output_path.status;
}

cairo_status_t
_cairo_path_stroke_to_path (const cairo_path_fixed_t		*path,
			    cairo_stroke_style_t *stroke_style,
			     cairo_path_fixed_t         *path_out,
			     cairo_matrix_t *ctm,
			     cairo_matrix_t *ctm_inverse,
			     double tolerance)
{
    path_ittr i;
    path_output_t output_path;

    cairo_status_t status = CAIRO_STATUS_SUCCESS;

    /* if we have an empty path, we're alredy done */
    if (cairo_path_head(path)->num_ops == 0)
	return CAIRO_STATUS_SUCCESS;

    if (stroke_style->num_dashes > 0) {
	if (path->has_curve_to) {
	    /* Dashing curved paths directly is tricky because I know of no easy way
	       to split a bezier curve at a particular length. So, for now, flatten the
	       path before dashing it */
	    cairo_path_fixed_t flat_path;
	    status = _cairo_path_fixed_init_flat_copy (&flat_path,
		    path,
		    tolerance);

	    if (unlikely(status)) {
		_cairo_path_fixed_fini (&flat_path);
		return status;
	    }

	    status = _cairo_path_dash_to_path (&flat_path, stroke_style, path_out, ctm, ctm_inverse);
	    _cairo_path_fixed_fini (&flat_path);
	    return status;
	} else {
	    return _cairo_path_dash_to_path (path, stroke_style, path_out, ctm, ctm_inverse);
	}
    }

    output_path.status = CAIRO_STATUS_SUCCESS;
    output_path.ctm = ctm;
    output_path.path = path_out;

    i.buf = i.first = cairo_path_head (path);
    i.op_index = 0;
    i.points = i.buf->points;

    while (!ittr_last (i)) {
	output_path.has_current_point = FALSE;
	_cairo_subpath_stroke_to_path (&i, stroke_style, &output_path, ctm_inverse);
    }

    return output_path.status;
}
