#include <stdio.h>
#include "tga.h"
#include "model.h"

void swap(int *a, int *b){
    int t = *a;
    *a = *b;
    *b = t;
}

int abs(int a){
    return (a >= 0) ? a : -a;
}

// Bresenham's line algorithm

void line (tgaImage *image, int x0, int y0, int x1, int y1, tgaColor color){
    int steep = 0;

    if (abs(y1 - y0) > abs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }

    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }

    int x;
    double y;
    double k = ((double)(y1 - y0))/(x1 - x0);

    for (x = x0, y = y0; x <= x1; ++x, y += k)
        if (steep)
            tgaSetPixel(image, y, x, color);
        else
            tgaSetPixel(image, x, y, color);
}

// Wireframe rendering of the model

void meshgrid(tgaImage *image, Model *model){
    int face, vert, i;
    int h = image->height;
    int w = image->width;
    tgaColor white = tgaRGB(255, 255, 255);
    Vec3 *p[3];
    for (face = 0; face < model->nface; ++face){
        for (vert = 0; vert < 3; ++vert){
            p[vert] = getVertex(model, face, vert);
        }
        for (i = 0; i < 3; ++i){
            line(image, ((*p[i])[0] + 1.0)*w/2, ((*p[i])[1] + 1.0)*h/2,
                ((*p[(i + 1)%3])[0] + 1.0)*w/2, ((*p[(i + 1)%3])[1] + 1.0)*h/2,
                white);
        }
    }
    tgaFlipVertically(image);
}

int main(int argc, char **argv){
    if (argc != 3){
        fprintf(stderr, "Usage: %s <objfile> <outfile>\n", argv[0]);
        return -1;
    }

    Model *model = NULL;
    model = loadFromObj(argv[1]);
    
    if (!model){
        perror("loadFromObj");
        freeModel(model);
        return -1;
    }

    tgaImage *image = NULL;
    image = tgaNewImage(800, 800, RGB);

    if (!image){
        perror("tgaNewImage");
        freeModel(model);
        tgaFreeImage(image);
        return -1;
    }

    meshgrid(image, model);

    if (-1 == tgaSaveToFile(image, argv[2])){
        perror("tgaSateToFile");
        freeModel(model);
        tgaFreeImage(image);
        return -1;
    }

    freeModel(model);
    tgaFreeImage(image);

    return 0;
}
