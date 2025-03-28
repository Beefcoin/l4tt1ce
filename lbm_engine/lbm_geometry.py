import numpy as np

""" This file contains some geometry helper methods to create geometry masks
In the Future i want to read geometry information from image data but atm  
this is a pretty good way to set up basic simulations for testing the implementation """

def create_triangle_mask(X, Y, center_x, center_y, base_width, height, direction='right'):
    mask = np.zeros_like(X, dtype=bool)

    if direction == 'right':
        x0 = center_x - base_width // 2
        x1 = center_x + height
        y0 = center_y - base_width // 2
        y1 = center_y + base_width // 2

        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                x, y = X[i, j], Y[i, j]
                if x0 <= x <= x1:
                    y_rel = (y - center_y) / (base_width / 2)
                    x_rel = (x - x0) / (x1 - x0)
                    if abs(y_rel) <= 1 - x_rel:
                        mask[i, j] = True

    elif direction == 'left':
        x0 = center_x + base_width // 2
        x1 = center_x - height
        y0 = center_y - base_width // 2
        y1 = center_y + base_width // 2

        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                x, y = X[i, j], Y[i, j]
                if x1 <= x <= x0:
                    y_rel = (y - center_y) / (base_width / 2)
                    x_rel = (x0 - x) / (x0 - x1)
                    if abs(y_rel) <= 1 - x_rel:
                        mask[i, j] = True

    return mask

def create_rectangle_mask(X, Y, center_x, center_y, width, height):
    mask = np.zeros_like(X, dtype=bool)

    x0 = center_x - width // 2
    x1 = center_x + width // 2
    y0 = center_y - height // 2
    y1 = center_y + height // 2

    mask[(X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)] = True

    return mask

def create_circle_mask(X, Y, center_x, center_y, radius):
    dist_squared = (X - center_x)**2 + (Y - center_y)**2
    mask = dist_squared <= radius**2
    return mask
