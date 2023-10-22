#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "tga.h"
#include "model.h"

#define IMAGE_WIDTH  800 // output image width
#define IMAGE_HEIGHT 800 // output image height

Vec3 light_dir = { -3. ,  1. , -1. }; // global light
Vec3 eye       = {  3. ,  1. ,  1. }; // camera position
Vec3 center    = {  0.6,  0. ,  0. }; // camera direction
Vec3 up        = {  0. ,  1. ,  0. }; // camera orientation

void swap(int *a, int *b){
    int t = *a;
    *a = *b;
    *b = t;
}

// Normalizing a vector

void normalize(Vec3 n){
    double norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; i++)
        n[i] /= norm;
}

// Matrix-Vector multiplication

void m2v(double m[4][4], double v[4], double result[4]){
    for (int i = 0; i < 4; i++){
        result[i] = 0.0;
        for (int j = 0; j < 4; j++)
            result[i] += m[i][j] * v[j];
    }
}

// Matrix-Matrix multiplication

void m2m(double m[4][4], double v[4][4], double result[4][4]){
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++){
            result[i][j] = 0;
            for (int k = 0; k < 4; k++)
                result[i][j] += m[i][k] * v[k][j];
        }
}

// Vectors-Vector multiplication

void cross(Vec3 a, Vec3 b, Vec3 c){
    c[0] = 0.; c[1] = 0.; c[2] = 0.;
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

// Computing the transition matrix "ModelView" for camera movement

void lookat(Vec3 eye, Vec3 center, Vec3 up, double ModelView[4][4]){
    Vec3 z, x, y;
    for (int i = 0; i < 3; i++)
        z[i] = eye[i] - center[i];
    normalize(z);

    cross(up, z, x); normalize(x);
    cross(z, x, y);  normalize(y);

    double Minv[4][4] = {
            {1.,0.,0.,0.},
            {0.,1.,0.,0.},
            {0.,0.,1.,0.},
            {0.,0.,0.,1.}};
    double Tr[4][4] = {
            {1.,0.,0.,0.},
            {0.,1.,0.,0.},
            {0.,0.,1.,0.},
            {0.,0.,0.,1.}};
    for (int i = 0; i < 3; i++){
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -center[i];
    }
    m2m(Minv, Tr, ModelView);
}

// Triangle drawing algorithm with z-buffer

void triangle(tgaImage *image, int sc[3][3], int uv[3][3], double intensity, int *zbuffer, Model *model){
    // Sort vertices by y coord
    if (sc[0][1] > sc[1][1])
        for (int k = 0; k < 3; k++){
            swap(&sc[0][k], &sc[1][k]);
            swap(&uv[0][k], &uv[1][k]);
        }
    if (sc[0][1] > sc[2][1])
        for (int k = 0; k < 3; k++){
            swap(&sc[0][k], &sc[2][k]);
            swap(&uv[0][k], &uv[2][k]);
        }
    if (sc[1][1] > sc[2][1])
        for (int k = 0; k < 3; k++){
            swap(&sc[1][k], &sc[2][k]);
            swap(&uv[1][k], &uv[2][k]);
        }
    int total_height = sc[2][1] - sc[0][1];
    for (int i = 0; i < total_height; i++){
        int second_half = i > sc[1][1] - sc[0][1] || sc[1][1] == sc[0][1];
        int segment_height = second_half ? sc[2][1] - sc[1][1] : sc[1][1] - sc[0][1];

        double alpha = (double) i / total_height;
        double beta  = (double)(i - (second_half ? sc[1][1] - sc[0][1] : 0)) / segment_height;

        int a[3], b[3], uva[3], uvb[3];
        for (int k = 0; k < 3; k++){
            a[k]   = rint(sc[0][k] + (sc[2][k] - sc[0][k]) * alpha);
            uva[k] = rint(uv[0][k] + (uv[2][k] - uv[0][k]) * alpha);
            
            b[k]   = rint(second_half ? sc[1][k] + (sc[2][k] - sc[1][k]) * beta
                                      : sc[0][k] + (sc[1][k] - sc[0][k]) * beta);
            uvb[k] = rint(second_half ? uv[1][k] + (uv[2][k] - uv[1][k]) * beta
                                      : uv[0][k] + (uv[1][k] - uv[0][k]) * beta);
        }
        if (a[0] > b[0])
            for (int k = 0; k < 3; k++){
                swap(  &a[k],   &b[k]);
                swap(&uva[k], &uvb[k]);
            }
        for (int j = a[0]; j <= b[0]; j++){
            double phi = (b[0] == a[0]) ? 1. : (double)(j - a[0]) / (double)(b[0] - a[0]);

            int p[3], uvp[3];
            for (int k = 0; k < 3; k++){
                p[k]   = rint(  a[k] + (  b[k] -   a[k]) * phi);
                uvp[k] = rint(uva[k] + (uvb[k] - uva[k]) * phi);
            }
            int idx = p[0] + p[1] * image->width;
            if (zbuffer[idx] < p[2]){
                zbuffer[idx] = p[2];

                Vec3 uv = {(double)(uvp[0]) / model->diffuse_map->width,
                           (double)(uvp[1]) / model->diffuse_map->height};

                tgaColor color = getDiffuseColor(model, &uv);
                color = tgaRGB(Red(color) * intensity, Green(color) * intensity, Blue(color) * intensity);

                tgaSetPixel(image, p[0], p[1], color);
            }
        }
    }
}

// Rendering of the model

void render(tgaImage *image, Model *model, int *zbuffer){
    int h = image->height, w = image->width, d = 255;

    double n[3];            // surface normal
    double v[4];            // object coords
    int sc[3][3];           // screen coords
    double wc[3][3];        // world coords;
    int uv[3][3];           // uv map;

    double  ModelView[4][4]; // transition matrix for camera movement
    lookat(eye, center, up, ModelView);

    double m[4] = {w/8., h/8., w*3./4., h*3./4.};
    double Viewport[4][4] = {
              {m[2]/2.,      0.,   0., m[0]+m[2]/2.},
              {     0., m[3]/2.,   0., m[1]+m[3]/2.},
              {     0.,      0., d/2.,         d/2.},
              {     0.,      0.,   0.,           1.}}; // transition matrix from world to camera
    
    double norm = sqrt(pow(eye[0] - center[0], 2) + pow(eye[1] - center[1], 2) + pow(eye[2] - center[2], 2));
    double Projection[4][4] = {{1., 0.,      0., 0.},
                               {0., 1.,      0., 0.},
                               {0., 0.,      1., 0.},
                               {0., 0.,-1./norm, 1.}}; // perspective projection
    
    Vec3 *p[3], *u[3];
    for (int face = 0; face < model->nface; ++face){
        for (int vert = 0; vert < 3; ++vert){
            p[vert] = getVertex(model, face, vert);
            u[vert] = getDiffuseUV(model, face, vert);
            v[0] = (*p[vert])[0];
            v[1] = (*p[vert])[1];
            v[2] = (*p[vert])[2];
            v[3] = 1.;
            
            // Viewport * Projection * View * Model * v.
            double r1[4], r2[4];
            m2v( ModelView,  v, r1);
            m2v(Projection, r1, r2);
            m2v(  Viewport, r2, r1);

            for (int k = 0; k < 3; k++){
                sc[vert][k] = r1[k]/r1[3];
                wc[vert][k] = (*p[vert])[k];
            }

            uv[vert][0] = (1. - (*u[vert])[0]) * model->diffuse_map->width;
            uv[vert][1] = (1. - (*u[vert])[1]) * model->diffuse_map->height; 
        }
        // surface normal calculation
        n[0] = (wc[2][1] - wc[0][1]) * (wc[1][2] - wc[0][2]) - (wc[2][2] - wc[0][2]) * (wc[1][1] - wc[0][1]);
        n[1] = (wc[2][2] - wc[0][2]) * (wc[1][0] - wc[0][0]) - (wc[2][0] - wc[0][0]) * (wc[1][2] - wc[0][2]); 
        n[2] = (wc[2][0] - wc[0][0]) * (wc[1][1] - wc[0][1]) - (wc[2][1] - wc[0][1]) * (wc[1][0] - wc[0][0]);

        normalize(n); normalize(light_dir);
        double intensity = n[0] * light_dir[0] + n[1] * light_dir[1] + n[2] * light_dir[2];

        if (intensity > 0)
            triangle(image, sc, uv, intensity, zbuffer, model);
    }
    tgaFlipVertically(image);    
}

int main(int argc, char **argv){
    int rm = 0;

    if (argc != 4){ fprintf(stderr, "Usage: %s <objfile> <diffusemap> <outfile>\n", argv[0]); return -1; }

    Model *model = NULL;
    model = loadFromObj(argv[1]);
    
    tgaImage *image = NULL;
    image = tgaNewImage(IMAGE_WIDTH, IMAGE_HEIGHT, RGB);

    int *zbuffer = NULL;
    zbuffer = malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(int));
    for (int i = 0; i < IMAGE_WIDTH * IMAGE_HEIGHT; i++)
        zbuffer[i] = INT_MIN;
    
    do{
        if (!model){
            perror("loadFromObj");
            rm = -1; break;
        }
        if (!image){
            perror("tgaNewImage");
            rm = -1; break;
        }
        if (!loadDiffuseMap(model, argv[2])){
            perror("loadDiffuseMap");
            rm = -1; break;
        }
        tgaFlipVertically(model->diffuse_map);
        tgaFlipHorizontally(model->diffuse_map);

        render(image, model, zbuffer);

        if (-1 == tgaSaveToFile(image, argv[3])){
            perror("tgaSateToFile");
            rm = -1; break;
        }
    } while(0);

    free(zbuffer);
    freeModel(model);
    tgaFreeImage(image);

    return rm;
}