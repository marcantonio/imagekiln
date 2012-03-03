
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <imagekiln.h>

// Fill color for coords outside of dest space.
#define FILL 0x00000000;

/*
 * Calculate a _reverse_ perspective transform matrix (perspective
 * transform, projective transform, homography, etc...).
 *
 *     Ax + By + C      Dx + Ev + F
 * u = -----------  v = -----------
 *     Gx + Hy + I      Gx + Hv + I
 *
 * Basically, a reverse mapping says: given the destination
 * coordinates (x,y) give me the cooresponding source coordinates.  A
 * forward mapping is also possible but there will be more "holes" in
 * the destination image and thus more information must be
 * interpolated.
 *
 * Homogeneous coordinate notation is used in these equations to
 * provide depth.  In Euclidean geometry there is no concept of
 * infinity, all 2D points are represented as (x,y).  In projective
 * geometry we need infinity to introduce depth.  Thus, a 2D point is
 * represented as a 3-vector (x',y',w).  Let x' = x/w where x is a
 * fixed value and x/0 is infinity.  The smaller w gets the larger the
 * value of x' becomes.  As w approaches 0, x' approaches infinity.
 * The homogeneous vector V = (x',y',w) = (x*w,y*w,w) represents the
 * actual point (x,y) = (x'/w,y'/w).  In practice, any non-zero number
 * can be used as w, but using w = 1 makes things simpler.
 *
 * For the mappings discussed, source point Ps = (u',v',q) and
 * destination point Pd = (x',y',w).  The forward projection matrix
 * will be shown as Msd and the reverse mapping shown as the adjoint
 * Mds.
 *
 * To transform any destination point to it's cooresponding source
 * point all you need to do is multiple the point by the inverse
 * mapping matrix:
 *
 * Ps = Pd * Mds
 *
 *                           /A D G\
 * (u',v',q)  =  (x',y',w) * |B E H|
 *                           \C F I/
 *
 * To solve the equation above we must solve 8 simultanious equations
 * for the coefficents used: A - H (I is always 1).  As descibed in
 * [Heckbert89] (pages 19-21) there is a quicker way when going from a
 * quadrilateral to a quadrilateral.
 *
 * First we compute the transform for a quadrilateral to a unit square
 * mapping (see Heckbert's paper for the simplified equations used
 * when we know it's a quad to square mapping).  Second, we compute
 * the adjoint (inverse) transform.  This adjoint transform is
 * actually a unit square to quadrilateral transform.  If we mulitply
 * the two matrices together (quad to square * square to quad) we get
 * a general quadilateral to quadrilateral transform matrix.
 *
 * This matrix becomes our Mds from above.  We then multiple each
 * homogenized destination point (homogenized == add a 1 on the end)
 * by the transform.  This gives us our homogeneous source point which
 * we must then normalize.
 *
 * The resulting point is then assigned a pixel color through linear
 * interpolation (area re-sampling).  This our anti-alias function
 * which gives the new image smoother lines and transitions.
 *
 * See Paul Heckbert's Master's thesis (sections 1 and 2):
 *   http://www.cs.cmu.edu/~ph/texfund/texfund.pdf
 *
 * Leptonica's affine transform page:
 *   http://www.leptonica.com/affine.html#RELATED-TRANSFORMS
 *
 * and this post on Stack Overflow:
 *   http://bit.ly/vxnpI5
 */
void perspective_transform(image_t *src_image, image_t *dst_image, vertex_t *v) {
	/*
	 * Four pairs of vertices for the transformation operation.  Points
	 * (u,v) denote the original quadrilateral.  Points (x,y) are the
	 * new locations for the transformed quadrilateral.
	 */
//	vertex_t vert[4];
//	vert[0].u = 148;
//	vert[0].v =  75;
//	vert[0].x =   0;
//	vert[0].y =   0;
//
//	vert[1].u = 615;
//	vert[1].v =  34;
//	vert[1].x = 639;
//	vert[1].y =   0;
//
//	vert[2].u = 585;
//	vert[2].v = 356;
//	vert[2].x = 639;
//	vert[2].y = 399;
//
//	vert[3].u =  70;
//	vert[3].v = 325;
//	vert[3].x =   0;
//	vert[3].y = 399;

	double_t xform[3][3];

	/* Calculate the transform needed. */
	quad_to_quad(v, xform);

	int x, y;
	pixel_t *ptr;
	pixel_t pval;
	double_t svec[3], dvec[3];

	/*
	 * For each point in dest space, find the cooresponding source space
	 * point and then find a pixel color to assign.
	 */
	for (y = 0; y < dst_image->height; y++) {
		ptr = dst_image->pixels + (y * dst_image->stride);
		for (x = 0; x < dst_image->width; x++) {
			// Make homogeneous.
			dvec[0] = x;
			dvec[1] = y;
			dvec[2] = 1.;

			// Do the transform.
			transform_point(dvec, xform, svec);

			// Normalize, i.e., (u,v) = (u'/q,v'/q)
			svec[0] /= svec[2];
			svec[1] /= svec[2];

			linear_interpolate(src_image->pixels, src_image->stride,
					dst_image->width, dst_image->height, svec[0], svec[1], &pval);

			*(ptr + x) = pval;
		}
	}
}

/*
 * Interpolate an color based on weighted area values.  Adapted from
 * Leptonica.
 */
void linear_interpolate(pixel_t *pixels, size_t stride, size_t width, size_t height,
		double_t x, double_t y, pixel_t *pval) {
	int32_t xpm, ypm, xp, yp, xf, yf;
	int32_t rval, gval, bval;
	uint32_t word00, word01, word10, word11;
	pixel_t *lines;
	pixel_t colorval = FILL;

	*pval = colorval;

	/* Skip if off the edge */
	if (x < 0.0 || y < 0.0 || x > width - 2.0 || y > height - 2.0)
		return;

	xpm = (int32_t)(16.0 * x + 0.5);
	ypm = (int32_t)(16.0 * y + 0.5);
	xp = xpm >> 4;
	yp = ypm >> 4;
	xf = xpm & 0x0f;
	yf = ypm & 0x0f;

	lines = pixels + (yp * stride);

	word00 = *(lines + xp);
	word10 = *(lines + xp + 1);
	word01 = *(lines + stride + xp);
	word11 = *(lines + stride + xp + 1);

#define L_RED_SHIFT  16
#define L_GREEN_SHIFT 8
#define L_BLUE_SHIFT  0

	rval = ((16 - xf) * (16 - yf) * ((word00 >> L_RED_SHIFT) & 0xff) +
			xf * (16 - yf) * ((word10 >> L_RED_SHIFT) & 0xff) +
			(16 - xf) * yf * ((word01 >> L_RED_SHIFT) & 0xff) +
			xf * yf * ((word11 >> L_RED_SHIFT) & 0xff) + 128) / 256;

	gval = ((16 - xf) * (16 - yf) * ((word00 >> L_GREEN_SHIFT) & 0xff) +
			xf * (16 - yf) * ((word10 >> L_GREEN_SHIFT) & 0xff) +
			(16 - xf) * yf * ((word01 >> L_GREEN_SHIFT) & 0xff) +
			xf * yf * ((word11 >> L_GREEN_SHIFT) & 0xff) + 128) / 256;

	bval = ((16 - xf) * (16 - yf) * ((word00 >> L_BLUE_SHIFT) & 0xff) +
			xf * (16 - yf) * ((word10 >> L_BLUE_SHIFT) & 0xff) +
			(16 - xf) * yf * ((word01 >> L_BLUE_SHIFT) & 0xff) +
			xf * yf * ((word11 >> L_BLUE_SHIFT) & 0xff) + 128) / 256;

	*pval = (rval << L_RED_SHIFT) | (gval << L_GREEN_SHIFT) | (bval << L_BLUE_SHIFT);
}

// Find determinant
#define DET(a,b,c,d) ((a)*(d)-(b)*(c))

/*
 * Compute a square to quad transform matrix for the four vertices
 * given.
 */
void square_to_quad(double_t p[][2], double_t sx[][3]) {
	double_t lx, ly;

	lx = p[0][0] - p[1][0] + p[2][0] - p[3][0];
	ly = p[0][1] - p[1][1] + p[2][1] - p[3][1];

	// If true, this is an affine mapping which is much quicker to
	// compute.
	if ((lx == 0) && (ly == 0)) {
		sx[0][0] = p[1][0] - p[0][0];
		sx[1][0] = p[2][0] - p[1][0];
		sx[2][0] = p[0][0];
		sx[0][1] = p[1][1] - p[0][1];
		sx[1][1] = p[2][1] - p[1][1];
		sx[2][1] = p[0][1];
		sx[0][2] = 0.;
		sx[1][2] = 0.;
		sx[2][2] = 1.;
		return;
	}

	double dx1, dx2, dy1, dy2, den;

	dx1 = p[1][0] - p[2][0];
	dx2 = p[3][0] - p[2][0];
	dy1 = p[1][1] - p[2][1];
	dy2 = p[3][1] - p[2][1];
	den = DET(dx1, dx2, dy1, dy2);

	double a, b, c, d, e, f, g, h;

	g = DET(lx,dx2,ly,dy2) / den;
	h = DET(dx1,lx,dy1,ly) / den;

	a = p[1][0] - p[0][0] + g * p[1][0];
	b = p[3][0] - p[0][0] + h * p[3][0];
	c = p[0][0];

	d = p[1][1] - p[0][1] + g * p[1][1];
	e = p[3][1] - p[0][1] + h * p[3][1];
	f = p[0][1];

	sx[0][0] = a; sx[0][1] = d; sx[0][2] = g;
	sx[1][0] = b; sx[1][1] = e; sx[1][2] = h;
	sx[2][0] = c; sx[2][1] = f; sx[2][2] = 1.;
}

/*
 * Using Heckbert's technique, find a quad to quad transform for the
 * vertices given.
 */
void quad_to_quad(vertex_t *v, double_t xform[][3]) {
	double_t md[3][3], dm[3][3], sm[3][3];
	double_t s_pts[4][2];
	double_t d_pts[4][2];
	int n;

	// Make the vertices easier to use.
	for (n = 0; n < 4; n++) {
		s_pts[n][0] = v[n].u;
		s_pts[n][1] = v[n].v;
		d_pts[n][0] = v[n].x;
		d_pts[n][1] = v[n].y;
	}

	// Find square to quad transform for the dst points.
	square_to_quad(d_pts, dm);

	// Invert for a quad to square transform.
	if (find_adjoint(dm, md) == 0.)
		// TODO:  Handle this
//		printf("Error: det == 0\n");
		;

	// Find square to quad transform for the src points.
	square_to_quad(s_pts, sm);

	// Combined them for the quad to quad transform.
	mmult(md, sm, xform);
}

/*
 * Finds a's adjoint matrix b.
 */
double find_adjoint(double a[][3], double b[][3]) {
	b[0][0] = DET(a[1][1], a[1][2], a[2][1], a[2][2]);
	b[1][0] = DET(a[1][2], a[1][0], a[2][2], a[2][0]);
	b[2][0] = DET(a[1][0], a[1][1], a[2][0], a[2][1]);
	b[0][1] = DET(a[2][1], a[2][2], a[0][1], a[0][2]);
	b[1][1] = DET(a[2][2], a[2][0], a[0][2], a[0][0]);
	b[2][1] = DET(a[2][0], a[2][1], a[0][0], a[0][1]);
	b[0][2] = DET(a[0][1], a[0][2], a[1][1], a[1][2]);
	b[1][2] = DET(a[0][2], a[0][0], a[1][2], a[1][0]);
	b[2][2] = DET(a[0][0], a[0][1], a[1][0], a[1][1]);
	return a[0][0]*b[0][0] + a[0][1]*b[0][1] + a[0][2]*b[0][2];
}

/*
 * Multiple two matrices.  (a * b = c)
 */
void mmult(double a[][3], double b[][3], double c[][3]) {
	c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
	c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
	c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];
	c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
	c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
	c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
	c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
	c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
	c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
}

/*
 * Do the actual point transform.  Homogeneous coords.
 */
void transform_point(double *dv, double h[][3], double *sv) {
	sv[0] = dv[0]*h[0][0] + dv[1]*h[1][0] + dv[2]*h[2][0];
	sv[1] = dv[0]*h[0][1] + dv[1]*h[1][1] + dv[2]*h[2][1];
	sv[2] = dv[0]*h[0][2] + dv[1]*h[1][2] + dv[2]*h[2][2];
}
