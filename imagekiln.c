
#include "imagekiln.h"

/* Convolution Matrix Functions */

// Gaussian blur
void blur(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride) {
  int kernel_blur[3][3] = { { 1, 2, 1 },
                            { 2, 4, 2 },
                            { 1, 2, 1 } };

  int factor = 16;
	int offset = 0;
  
  convolve(src_pixels, dst_pixels, width, height, stride, kernel_blur, factor, offset);
}

// Edge detect
void edge_detect(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride) {
  int kernel_sobel1[3][3] = { { -1, -2, -1 },
                              {  0,  0,  0 },
                              {  1,  2,  1 } };
  
  int kernel_sobel2[3][3] = { { -1,  0,  1 },
                              { -2,  0,  2 },
                              { -1,  0,  1 } };
  
  int factor = 1;
  int offset = 0;
  int threshold = 80;
  
  convolve_multi_kernel_gray(src_pixels, dst_pixels, width, height, stride, 
                             kernel_sobel1, kernel_sobel2, factor, offset, threshold);
}

/* Histogram functions */

void histogram_equalize_color_to_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride) {
    image_copy(src_pixels, dst_pixels, height, stride);
    
    histogram_equalize_channel(src_pixels, dst_pixels, width, height, stride, GRAY);
}

void histogram_equalize_rgb(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride) {
    image_copy(src_pixels, dst_pixels, height, stride);
    histogram_equalize_channel(src_pixels, dst_pixels, width, height, stride, RED);
    histogram_equalize_channel(src_pixels, dst_pixels, width, height, stride, GREEN);
    histogram_equalize_channel(src_pixels, dst_pixels, width, height, stride, BLUE);    
}

/*
 * Assumption: 2^32 is enough to store the sum of the occurrences of each gray
 * level.
 *
 * Note: wikipedia indicates subtracting MINIMUM(CDF) from area and also from
 * current CDF value, when computing, but due to the fact that an implementation
 * I saw
 * (http://code.google.com/p/simple-iphone-image-processing/source/browse/trunk/Classes/Image.mm)
 * and A book (http://homepages.inf.ed.ac.uk/rbf/BOOKS/PHILLIPS/) seem to ignore
 * that and the reason for it isn't clear in the wikipedia article, so I'm
 * ignoring it here.
 */

void histogram_equalize_channel(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride, channels channel) {
  
    /*
     * Basic Algorithm
     * 1. Loop through source image.
     *  a. calculate channel value for each pixel
     *  b. Update histogram for each gray value.
     * 2. Loop through histogram
     *  a. calculate and store the sum of the histogram at each channel level.
     * 3. loop through source image again
     *  a. using the CDF and histogram computed in the first loop, apply the
     *     transform for the channel value of this pixel to each color pixel and create the
     *     destination image.
     */

    // Calculate histogram
    uint32_t histogram[UINT8_MAX + 1] = { 0 };
    create_image_histogram(src_pixels, width, height, stride, histogram, channel);

    // Calculate CDF
    uint32_t cdf[UINT8_MAX + 1] = { 0 };
    cdf[0] = histogram[0];
    for (int value_level = 1; value_level < UINT8_MAX + 1; value_level++) {
        cdf[value_level] = cdf[value_level - 1] + histogram[value_level];
    }

    // Convert on selected channel using calculated conversion factor.
    equalize_histogram_on_channel(src_pixels, dst_pixels, width, height, stride, cdf, channel);    
}

void create_image_histogram(void *src_pixels, size_t width, size_t height, size_t stride, uint32_t *histogram, channels channel) {

    value_function_ptr value_function = get_value_function_pointer(channel);
    
    for (int h = 0; h < height - 1; h++) {
        uint32_t *src_row = src_pixels + (stride * h);
        for (int w = 0; w < width - 1; w++) {
            uint8_t value = value_function(src_row[w]);
            histogram[value] += 1;
        }
    }
    
}

void equalize_histogram_on_channel(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride, uint32_t *cdf, channels channel) {

    double_t conversion_constant = get_conversion_constant(width, height);
    
    value_function_ptr value_function = get_value_function_pointer(channel);
    
    pixel_convert_function_ptr pixel_convert_function = get_pixel_convert_function_pointer(channel);

    for (int h = 0; h < height - 1; h++) {
        uint32_t *src_row = src_pixels + (stride * h);
        uint32_t *dst_row = dst_pixels + (stride * h);
        for (int w = 0; w < width - 1; w++) {
            uint32_t old_value = value_function(src_row[w]);
            uint32_t new_value = conversion_constant * cdf[old_value];
            dst_row[w] = pixel_convert_function(dst_row[w], old_value, new_value);
        }
    }
}

double_t get_conversion_constant(size_t width, size_t height) {
    // Calculate constant for conversion based on desired number of value levels and the area of the image
    uint32_t image_area = width * height;
    uint32_t values_in_output = UINT8_MAX + 1;
    return (double_t)values_in_output/(double_t)image_area;
}

/* Util functions */
void image_copy(void *src_pixels, void *dst_pixels, size_t height, size_t stride) {
    size_t dataSize = stride * height;
    memcpy(dst_pixels, src_pixels, dataSize);
}

value_function_ptr get_value_function_pointer(channels channel) {
    switch (channel) {
        case RED:
            return red_value;
        case GREEN:
            return green_value;
        case BLUE:
            return blue_value;
        case ALPHA:
            return alpha_value;
        case GRAY:
            return gray_value;
        default:
            return NULL;
    }
}

pixel_convert_function_ptr get_pixel_convert_function_pointer(channels channel) {
    switch (channel) {
        case RED:
            return pixel_convert_red;
        case GREEN:
            return pixel_convert_green;
        case BLUE:
            return pixel_convert_blue;
        case ALPHA:
            return pixel_convert_alpha;
        case GRAY:
            return pixel_convert_gray;
        default:
            return NULL;
    }
}

/*
 * Assumes 32 bit unsigned integer passed representing RGBA 8 bits each.  Multiples 
 * each value by passed factor.
 */
uint32_t multiply_channels_with_value(uint32_t rgba_pixel, double_t conversion_value) {
    uint8_t red, green, blue, alpha = 0;
    red = red_value(rgba_pixel) * conversion_value;
    green = green_value(rgba_pixel) * conversion_value;
    blue = blue_value(rgba_pixel) * conversion_value;
    alpha = alpha_value(rgba_pixel);
    return pixel_from_components(red, green, blue, alpha);
}

/* 
 * Assumes a 32 bit unsigned integer passed which represents RGBA of a pixel (8bits each)
 * Converts by using simple method of:
 * 30% R, 59% G, 11% B and I ignore anything beyond 8 bits of Blue. 
 */
uint8_t gray_value(uint32_t rgb_pixel) {
    uint8_t red_component, green_component, blue_component = 0;

    red_component = red_value(rgb_pixel) * 0.30f;
    green_component = green_value(rgb_pixel) * 0.59f;
    blue_component = blue_value(rgb_pixel) * 0.11f;

    return red_component + green_component + blue_component;
}

uint8_t red_value(uint32_t rgb_pixel) {
    return (uint8_t)(rgb_pixel & 0x000000ff);
}

uint8_t green_value(uint32_t rgb_pixel) {
    return (uint8_t)((rgb_pixel & 0x0000ff00) >> 8);
}

uint8_t blue_value(uint32_t rgb_pixel) {
    return (uint8_t)((rgb_pixel & 0x00ff0000) >> 16);
}

uint8_t alpha_value(uint32_t rgb_pixel) {
    return (uint8_t)((rgb_pixel & 0xff000000) >> 24);
}

uint32_t pixel_from_components(uint8_t red, uint8_t green, uint8_t blue, uint8_t alpha) {
    return ((alpha << 24) | (blue << 16) | (green << 8) | red);
}

uint32_t pixel_convert_red(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value) {
    return pixel_from_components(new_value, green_value(rgba_pixel), blue_value(rgba_pixel), alpha_value(rgba_pixel));
}

uint32_t pixel_convert_green(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value) {
    return pixel_from_components(red_value(rgba_pixel), new_value, blue_value(rgba_pixel), alpha_value(rgba_pixel));
}
uint32_t pixel_convert_blue(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value) {
    return pixel_from_components(red_value(rgba_pixel), green_value(rgba_pixel), new_value, alpha_value(rgba_pixel));
}
uint32_t pixel_convert_alpha(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value) {
    return pixel_from_components(red_value(rgba_pixel), green_value(rgba_pixel), blue_value(rgba_pixel), new_value);
}
uint32_t pixel_convert_gray(uint32_t rgba_pixel, uint8_t old_value, uint8_t new_value) {
    double_t conversion_value = (double_t)new_value/(double_t)old_value;
    return multiply_channels_with_value(rgba_pixel, conversion_value);
}

void grayscale_image(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride) {
	uint8_t r, g, b;
	int gray_value;
  
	for (int h = 0; h < height; h++) {
		uint32_t *src_row = (uint32_t *)src_pixels;
		uint32_t *dst_row = (uint32_t *)dst_pixels;
    
		int w;
		for (w = 0; w < width; w++) {
      //			a = src_row[w] >> 24;
			b = (src_row[w] & 0x00ff0000) >> 16;
			g = (src_row[w] & 0x0000ff00) >> 8;
			r = src_row[w] & 0x000000ff;
      
			gray_value = (r * GS_RED_FACTOR) + (g * GS_GREEN_FACTOR) + (b * GS_BLUE_FACTOR);
      //			dst_row[w] = (a << 24) | (gray_value << 16) | (gray_value << 8) | gray_value;
			dst_row[w] = (0 << 24) | (gray_value << 16) | (gray_value << 8) | gray_value;
		}
		src_pixels = (char *) src_pixels + stride;
		dst_pixels = (char *) dst_pixels + stride;
	}
}
