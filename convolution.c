#include <stdlib.h>
#include <android/log.h>

#define  LOG_TAG    "libnativekiln"
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)

#include "imagekiln.h"

void convolve(void *src_pixels, void *dst_pixels, size_t width, size_t height, 
              size_t stride, int kernel[][3], int factor, int offset) {
	uint32_t *src_row, *dst_row, *prev_row, *next_row;
	int16_t sR, sG, sB, sA;
  
	int h;
	for (h = 1; h < height - 1; h++) {
		src_pixels = (char *)src_pixels + stride;
		dst_pixels = (char *)dst_pixels + stride;
    
		src_row = (uint32_t *)src_pixels;
		dst_row = (uint32_t *)dst_pixels;
    
		prev_row = (uint32_t *)((char *)src_row - stride);
		next_row = (uint32_t *)((char *)src_row + stride);
    
		// Each pixel
		int w;
		for (w = 1; w < width - 1; w++) {
			int i, x;
			sR = sG = sB = 0;
			sA = src_row[w] >> 24;
      
			// Image is stored upside down.
			for (i = w - 1, x = 0; i <= w + 1; i++, x++) {
        // bottom row
				sB += ((prev_row[i] & 0x00ff0000) >> 16) * kernel[2][x];
				sG += ((prev_row[i] & 0x0000ff00) >> 8) * kernel[2][x];
				sR += (prev_row[i] & 0x000000ff) * kernel[2][x];
        
        // center row
				sB += ((src_row[i] & 0x00ff0000) >> 16) * kernel[1][x];
				sG += ((src_row[i] & 0x0000ff00) >> 8) * kernel[1][x];
				sR += (src_row[i] & 0x000000ff) * kernel[1][x];
        
        // top row
				sB += ((next_row[i] & 0x00ff0000) >> 16) * kernel[0][x];
				sG += ((next_row[i] & 0x0000ff00) >> 8) * kernel[0][x];
				sR += (next_row[i] & 0x000000ff) * kernel[0][x];
			}
      
			sR = abs((float)sR / factor) + offset;
			sG = abs((float)sG / factor) + offset;
			sB = abs((float)sB / factor) + offset;
      
      // clamp
			if (sR > 255) {
				sR = 255;
			} else if (sR < 0) {
				sR = 0;
			}
      
			if (sG > 255) {
				sG = 255;
			} else if (sG < 0) {
				sG = 0;
			}
      
			if (sB > 255) {
				sB = 255;
			} else if (sB < 0) {
				sB = 0;
			}
      
      //			dst_row[w] = ((sA << 24) | (sB << 16) | (sG << 8) | sR);
			dst_row[w] = ((0 << 24) | (sB << 16) | (sG << 8) | sR);
		}
	}
}

void convolve_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, 
                   size_t stride, int kernel[][3], int factor, int offset) {
	uint32_t *src_row, *dst_row, *prev_row, *next_row;
	int16_t pixel_value;
      
	int h;
	for (h = 1; h < height - 1; h++) {
		src_pixels = (char *)src_pixels + stride;
		dst_pixels = (char *)dst_pixels + stride;
    
		src_row = (uint32_t *)src_pixels;
		dst_row = (uint32_t *)dst_pixels;
    
		prev_row = (uint32_t *)((char *)src_row - stride);
		next_row = (uint32_t *)((char *)src_row + stride);
    
		// Each pixel
		int w;
		for (w = 1; w < width - 1; w++) {
			int i, x;
			pixel_value = 0;
      
			// Image is stored upside down.
			for (i = w - 1, x = 0; i <= w + 1; i++, x++) {
        // bottom row
				pixel_value += (prev_row[i] & 0x000000ff) * kernel[2][x];
        
        // center row
				pixel_value += (src_row[i] & 0x000000ff) * kernel[1][x];
        
        // top row
				pixel_value += (next_row[i] & 0x000000ff) * kernel[0][x];
			}
      
      pixel_value = abs((float)pixel_value / factor) + offset;
      
      // clamp; must be signed!!!!
			if (pixel_value > 255) {
				pixel_value = 255;
			} else if ((int16_t)pixel_value < 0) {
				pixel_value = 0;
			}

			dst_row[w] = ((0 << 24) | (pixel_value << 16) | (pixel_value << 8) | pixel_value);
		}
	}
}

void convolve_multi_kernel_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride, 
                                int kernel1[][3], int kernel2[][3], int factor, int offset, int threshold) {
	uint32_t *src_row, *dst_row, *prev_row, *next_row;
  //	uint16_t pixel_value1, pixel_value2;
	int16_t pixel_value1, pixel_value2;
  
	int h;
	for (h = 1; h < height - 1; h++) {
		src_pixels += stride;
		dst_pixels += stride;
    
		src_row = (uint32_t *)src_pixels;
		dst_row = (uint32_t *)dst_pixels;
    
		prev_row = (uint32_t *)((unsigned char *)src_row - stride);
		next_row = (uint32_t *)((unsigned char *)src_row + stride);
    
		// Each pixel
		int w;
		for (w = 1; w < width - 1; w++) {
			int i, x;
			pixel_value1 = 0;
			pixel_value2 = 0;
      
			// Image is stored upside down.
			for (i = w - 1, x = 0; i <= w + 1; i++, x++) {
        // bottom row
				pixel_value1 += (prev_row[i] & 0x000000ff) * kernel1[2][x];
				pixel_value2 += (prev_row[i] & 0x000000ff) * kernel2[2][x];
        
        // center row
				pixel_value1 += (src_row[i] & 0x000000ff) * kernel1[1][x];
				pixel_value2 += (src_row[i] & 0x000000ff) * kernel2[1][x];
        
        // top row
				pixel_value1 += (next_row[i] & 0x000000ff) * kernel1[0][x];
				pixel_value2 += (next_row[i] & 0x000000ff) * kernel2[0][x];
			}
      
      pixel_value1 = abs((float)pixel_value1 / factor) + offset;
      pixel_value2 = abs((float)pixel_value2 / factor) + offset;
      
      int16_t pixel_value = sqrt((pixel_value1*pixel_value1) + (pixel_value2*pixel_value2));
      
      pixel_value = (pixel_value > threshold) ? 255 : 0;
      
			dst_row[w] = ((0 << 24) | (pixel_value << 16) | (pixel_value << 8) | pixel_value);
		}
	}
}
