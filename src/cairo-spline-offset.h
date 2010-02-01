typedef cairo_point_double_t point_t;
typedef struct _knots {
	    point_t a,b,c,d;
} knots_t;

void print_knot(knots_t k);

struct curve_list {
	knots_t curve;
	double dist;
	struct curve_list *next;
	struct curve_list *prev;
};

typedef struct curve_list curve_list_t;

void
curve_offset(knots_t self, double offset, double tolerance, void (*curve_to_fn)(void *path, knots_t k), void *closure);
