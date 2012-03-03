
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __ANDROID__
	typedef double double_t;
#endif

/* Types */
enum channel_names {
    RED = 0,
    GREEN = 1,
    BLUE = 2,
    ALPHA = 3,
    GRAY = 4
};

typedef enum channel_names channels;

// pointer to function that takes pixel and returns a single 8 bit value
// Used to generalize stripping out red/gree/blue/alpha or converting to gray.
typedef uint8_t (*value_function_ptr)(uint32_t rgba_pixel);

// pointer to a function that takes pixel, old value, new value and returns a new pixel
// used to generalize converting values based on channel.
typedef uint32_t (*pixel_convert_function_ptr)(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value);

typedef uint32_t pixel_t;

// Grayscale constants
#define GS_RED_FACTOR   0.299f
#define GS_GREEN_FACTOR 0.587f
#define GS_BLUE_FACTOR  0.114f

// Convenient type to pass around important image data.
typedef struct {
  uint32_t width;
  uint32_t height;
  uint32_t size;
  uint32_t stride;
  pixel_t *pixels;
} image_t;

/*
 * Represents the source coordinates (u,v) and the destination
 * coordinates (x,y).
 */
typedef struct {
  double u, v;
  double x, y;
} vertex_t;

/* prototypes */
void image_copy(void *src_pixels, void *dst_pixels, size_t height, size_t stride);

uint8_t red_value(uint32_t rgb_pixel);
uint8_t green_value(uint32_t rgb_pixel);
uint8_t blue_value(uint32_t rgb_pixel);
uint8_t alpha_value(uint32_t rgb_pixel);

uint8_t gray_value(uint32_t rgb_pixel);

value_function_ptr get_value_function_pointer(channels channel);

uint32_t pixel_from_components(uint8_t red, uint8_t green, uint8_t blue, uint8_t alpha);

uint32_t pixel_convert_red(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value);
uint32_t pixel_convert_green(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value);
uint32_t pixel_convert_blue(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value);
uint32_t pixel_convert_alpha(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value);
uint32_t pixel_convert_gray(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value);

pixel_convert_function_ptr get_pixel_convert_function_pointer(channels channel);

uint32_t multiply_channels_with_value(uint32_t rgba_pixel, double_t conversion_value);

/* histogram prototypes */
double_t get_conversion_constant(size_t width, size_t height);
void create_image_histogram(void *src_pixels, size_t width, size_t height, size_t stride, uint32_t *histogram, channels channel);
void equalize_histogram_on_channel(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride, uint32_t *cdf, channels channel);
void histogram_equalize_channel(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride, channels channel);

/* Perspective tranform prototypes */
void quad_to_quad(vertex_t *v, double_t xform[][3]);
void square_to_quad(double_t p[][2], double_t sx[][3]);
double_t find_adjoint(double_t a[][3], double_t b[][3]);
void mmult(double_t a[][3], double_t b[][3], double_t c[][3]);
void transform_point(double_t *dv, double_t h[][3], double_t *sv);
void linear_interpolate(pixel_t *pixels, size_t stride, size_t width, size_t height, double_t x, double_t y, pixel_t *pval);

/* Generic convolution function. Kernel defined in the calling function. */
void convolve(void *src_pixels, void *dst_pixels, size_t width, size_t height, 
              size_t stride, int kernel[][3], int factor, int offset);
/* Convole a grayscale image */
void convolve_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, 
                   size_t stride, int kernel[][3], int factor, int offset);
/* Convole a grayscale image with two kernels.  The two results are simply added together.*/
void convolve_multi_kernel_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride, 
                                int kernel1[][3], int kernel2[][3], int factor, int offset, int threshold);

/* prototypes and types a user of the library would need to know and use */
void blur(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void histogram_equalize_color_to_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void histogram_equalize_rgb(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void perspective_transform(image_t *src_image, image_t *dst_image, vertex_t *v);
void grayscale_image(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void edge_detect(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
