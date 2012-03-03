
/* prototypes and types a user of the library would need to know and use */

void blur(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void histogram_equalize_color_to_gray(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void histogram_equalize_rgb(void *src_pixels, void *dst_pixels, size_t width, size_t height, size_t stride);
void perspective_transform(image_t *src_image, image_t *dst_image, vertex_t *v);
