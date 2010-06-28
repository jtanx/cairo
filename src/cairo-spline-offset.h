#ifndef CAIRO_SPLINE_OFFSET_H
#define CAIRO_SPLINE_OFFSET_H

cairo_status_t
_cairo_spline_offset (const cairo_spline_knots_double_t *k, double offset, double tolerance, cairo_status_t (*curve_fn)(void *, const cairo_spline_knots_double_t *), void *closure);

#endif /*CAIRO_DIRECTFB_H*/

