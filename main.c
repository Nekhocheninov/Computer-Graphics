#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "tga.h"
#include "model.h"

void swap(int *a, int *b){
    int t = *a;
    *a = *b;
    *b = t;
}

int *zbuffer = NULL;

// Triangle drawing algorithm with z-buffer

void triangle(tgaImage *image, int x0, int y0, int z0,
                               int x1, int y1, int z1,
                               int x2, int y2, int z2, tgaColor color, int *zbuffer){
    // Sort vertices by y coord
    if (y0 > y1){
        swap(&x0, &x1);
        swap(&y0, &y1);
        swap(&z0, &z1);
    }
    if (y0 > y2){
        swap(&x0, &x2);
        swap(&y0, &y2);
        swap(&z0, &z2);
    }
    if (y1 > y2){
        swap(&x1, &x2);
        swap(&y1, &y2);
        swap(&z1, &z2);
    }

    int total_height = y2 - y0;
    for (int i = 0; i < total_height; i++){
        int second_half = i > y1 - y0 || y1 == y0;
        int segment_height = second_half ? y2 - y1 : y1 - y0;
        double alpha = (double) i / total_height;
        double beta  = (double)(i - (second_half ? y1 - y0 : 0)) / segment_height; // be careful: with above conditions no division by zero here

        int ax = rint(x0 + (x2 - x0) * alpha);
        int ay = rint(y0 + (y2 - y0) * alpha);
        int az = rint(z0 + (z2 - z0) * alpha);
        int bx = rint(second_half ? x1 + (x2 - x1)*beta : x0 + (x1 - x0) * beta);
        int by = rint(second_half ? y1 + (y2 - y1)*beta : y0 + (y1 - y0) * beta);
        int bz = rint(second_half ? z1 + (z2 - z1)*beta : z0 + (z1 - z0) * beta);
        if (ax > bx){
            swap(&ax, &bx);
            swap(&ay, &by);
            swap(&az, &bz);
        }

        for (int j = ax; j <= bx; j++){
            double phi = (bx == ax) ? 1. : (double)(j - ax) / (double)(bx - ax);
            int px = rint(ax + (bx - ax) * phi);
            int py = rint(ay + (by - ay) * phi);
            int pz = rint(az + (bz - az) * phi);
            int idx = px + py * image->width;
            if (zbuffer[idx] < pz){
                zbuffer[idx] = pz;
                tgaSetPixel(image, px, py, color);
            }
        }
    }
}

// Rendering of the model

void render(tgaImage *image, Model *model){
    int face, vert;
    int h = image->height;
    int w = image->width;
    int d = 255; // depth for z buffer

    zbuffer = malloc(h * w * sizeof(int));
    for (int i = 0; i < h * w; i++)
        zbuffer[i] = INT_MIN;

    Vec3 *p[3];
    double n[3];
    int sc[3][3]; // screen coords
    double wc[3][3]; // world coords;

    for (face = 0; face < model->nface; ++face){
        for (vert = 0; vert < 3; ++vert){
            p[vert] = getVertex(model, face, vert);
            sc[vert][0] = ((*p[vert])[0] + 1.0)*w/2;
            sc[vert][1] = ((*p[vert])[1] + 1.0)*h/2;
            sc[vert][2] = ((*p[vert])[2] + 1.0)*d/2;
            wc[vert][0] = (*p[vert])[0];
            wc[vert][1] = (*p[vert])[1];
            wc[vert][2] = (*p[vert])[2];

        }
        n[0] = (wc[2][1] - wc[0][1]) * (wc[1][2] - wc[0][2]) - (wc[2][2] - wc[0][2]) * (wc[1][1] - wc[0][1]);
        n[1] = (wc[2][2] - wc[0][2]) * (wc[1][0] - wc[0][0]) - (wc[2][0] - wc[0][0]) * (wc[1][2] - wc[0][2]); 
        n[2] = (wc[2][0] - wc[0][0]) * (wc[1][1] - wc[0][1]) - (wc[2][1] - wc[0][1]) * (wc[1][0] - wc[0][0]);

        double norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

        n[0] = n[0] / norm;
        n[1] = n[1] / norm;
        n[2] = n[2] / norm;

        double intensity = n[0] * 0. + n[1] * 0. + n[2] * (-1.);
        
        if (intensity > 0)
            triangle(image, sc[0][0], sc[0][1], sc[0][2],
                            sc[1][0], sc[1][1], sc[1][2],
                            sc[2][0], sc[2][1], sc[2][2],
                            tgaRGB(intensity*255, intensity*255, intensity*255), zbuffer);
    }
    free(zbuffer);
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
    
    render(image, model);

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
