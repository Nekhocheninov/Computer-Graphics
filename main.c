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

void normalize(double n[3]){
    double norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    for (int i = 0; i < 3; i++)
        n[i] /= norm;
}

void m2v(double m[4][4], double v[4], double result[4]){
    for (int i = 0; i < 4; i++){
        result[i] = 0.0;
        for (int j = 0; j < 4; j++)
            result[i] += m[i][j] * v[j];
    }
}

void m2m(double m[4][4], double v[4][4], double result[4][4]){
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            result[i][j] = 0;
            for (int k = 0; k < 4; k++)
                result[i][j] += m[i][k] * v[k][j];
        }
}

void cross(Vec3 a, Vec3 b, Vec3 c){
    c[0] = 0.; c[1] = 0.; c[2] = 0.;
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

int *zbuffer = NULL;

Vec3 light_dir = {1., 0., 1.};
Vec3 eye = {-3.5, -1., -5.};
Vec3 center = {0.6, 0., 0.};
Vec3 up = {0., 1., 0.};

double ModelView[4][4];

void lookat(Vec3 eye, Vec3 center, Vec3 up) {
    Vec3 z, x, y;
    for (int i = 0; i < 3; i++)
        z[i] = eye[i] - center[i];
    normalize(z);

    cross(up, z, x);
    normalize(x);

    cross(z, x, y);
    normalize(y);
    double Minv[4][4] = {{1., 0., 0., 0.},
                         {0., 1., 0., 0.},
                         {0., 0., 1., 0.},
                         {0., 0., 0., 1.}};
    double Tr[4][4] =   {{1., 0., 0., 0.},
                         {0., 1., 0., 0.},
                         {0., 0., 1., 0.},
                         {0., 0., 0., 1.}};
    for (int i = 0; i < 3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -center[i];
    }
    m2m(Minv, Tr, ModelView);
}

// Triangle drawing algorithm with z-buffer

void triangle(tgaImage *image, int x0, int y0, int z0,
                               int x1, int y1, int z1,
                               int x2, int y2, int z2,
                               int uvx0, int uvy0,
                               int uvx1, int uvy1,
                               int uvx2, int uvy2,
                               double intensity, int *zbuffer, Model *model){
    // Sort vertices by y coord
    if (y0 > y1){
        swap(&x0, &x1);
        swap(&y0, &y1);
        swap(&z0, &z1);
        swap(&uvx0, &uvx1);
        swap(&uvy0, &uvy1);
    }
    if (y0 > y2){
        swap(&x0, &x2);
        swap(&y0, &y2);
        swap(&z0, &z2);
        swap(&uvx0, &uvx2);
        swap(&uvy0, &uvy2);
    }
    if (y1 > y2){
        swap(&x1, &x2);
        swap(&y1, &y2);
        swap(&z1, &z2);
        swap(&uvx1, &uvx2);
        swap(&uvy1, &uvy2);
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
        int uvax = rint(uvx0 + (uvx2 - uvx0)*alpha);
        int uvbx = rint(second_half ? uvx1 + (uvx2-uvx1)*beta : uvx0 + (uvx1-uvx0)*beta);
        int uvay = rint(uvy0 + (uvy2 - uvy0)*alpha);
        int uvby = rint(second_half ? uvy1 + (uvy2-uvy1)*beta : uvy0 + (uvy1-uvy0)*beta);
        if (ax > bx){
            swap(&ax, &bx);
            swap(&ay, &by);
            swap(&az, &bz);

            swap(&uvax, &uvbx);
            swap(&uvay, &uvby);
        }

        for (int j = ax; j <= bx; j++){
            double phi = (bx == ax) ? 1. : (double)(j - ax) / (double)(bx - ax);
            int px = rint(ax + (bx - ax) * phi);
            int py = rint(ay + (by - ay) * phi);
            int pz = rint(az + (bz - az) * phi);
            int idx = px + py * image->width;
            int uvpx = rint(abs(uvax + (uvbx - uvax) * phi));
            int uvpy = rint(abs(uvay + (uvby - uvay) * phi));
            if (zbuffer[idx] < pz){
                zbuffer[idx] = pz;
                Vec3 temp = {0.0, 0.0, 0.0};
                Vec3 *uv = &temp;
                double width = model->diffuse_map->width;
                double height = model->diffuse_map->height;
                (*uv)[0] = uvpx/width;
                (*uv)[1] = uvpy/height;

                tgaColor color = getDiffuseColor(model, uv);

                color = tgaRGB(Red(color) * intensity, Green(color) * intensity, Blue(color) * intensity);
                tgaSetPixel(image, px, py, color); // tgaSetPixel(image, j, u0+i, color);
            }
        }
    }
}

// Rendering of the model

void render(tgaImage *image, Model *model){
    int face, vert;
    int h = image->height;
    int w = image->width;
    int d = 255;

    zbuffer = malloc(h * w * sizeof(int));
    for (int i = 0; i < h * w; i++)
        zbuffer[i] = INT_MIN;

    Vec3 *p[3];
    Vec3 *uv[3];
    double n[3];
    int sc[3][3]; // screen coords
    double wc[3][3]; // world coords;
    int tc[3][2]; // uv map;

    lookat(eye, center, up);
    
    double ViewPort[4] = {w/8., h/8., w*3./4., h*3./4.};
    double m[4][4] =  {{ViewPort[2]/2., 0., 0., ViewPort[0]+ViewPort[2]/2.},
                        {0., ViewPort[3]/2., 0., ViewPort[1]+ViewPort[3]/2.},
                        {0., 0., d/2., d/2.},
                        {0., 0., 0., 1.}};
    Vec3 camera;
    for (int i = 0; i < 3; i++)
        camera[i] = (eye[i] - center[i]);
    double norm = sqrt(camera[0]*camera[0] + camera[1]*camera[1] + camera[2]*camera[2]);
    double pj[4][4] = {{1., 0., 0., 0.},
                       {0., 1., 0., 0.},
                       {0., 0., 1., 0.},
                       {0., 0.,-1./norm, 1.}};
    double v[4];

    for (face = 0; face < model->nface; ++face){
        for (vert = 0; vert < 3; ++vert){
            p[vert] = getVertex(model, face, vert);
            uv[vert] = getDiffuseUV(model, face, vert);
            v[0] = (*p[vert])[0];
            v[1] = (*p[vert])[1];
            v[2] = (*p[vert])[2];
            v[3] = 1.;
            
            double r1[4], r2[4];
            m2v(ModelView, v, r1);
            m2v(pj, r1, r2);
            m2v(m, r2, r1);

            
            sc[vert][0] = r1[0]/r1[3];
            sc[vert][1] = r1[1]/r1[3];
            sc[vert][2] = r1[2]/r1[3];

            wc[vert][0] = (*p[vert])[0];
            wc[vert][1] = (*p[vert])[1];
            wc[vert][2] = (*p[vert])[2];
            
            tc[vert][0] = (1. - (*uv[vert])[0]) * model->diffuse_map->width;
            tc[vert][1] = (1. - (*uv[vert])[1]) * model->diffuse_map->height; 
        }
        n[0] = (wc[2][1] - wc[0][1]) * (wc[1][2] - wc[0][2]) - (wc[2][2] - wc[0][2]) * (wc[1][1] - wc[0][1]);
        n[1] = (wc[2][2] - wc[0][2]) * (wc[1][0] - wc[0][0]) - (wc[2][0] - wc[0][0]) * (wc[1][2] - wc[0][2]); 
        n[2] = (wc[2][0] - wc[0][0]) * (wc[1][1] - wc[0][1]) - (wc[2][1] - wc[0][1]) * (wc[1][0] - wc[0][0]);

        normalize(n); normalize(light_dir);
    
        double intensity = n[0] * light_dir[0] + n[1] * light_dir[1] + n[2] * light_dir[2];

        if (intensity > 0)
            triangle(image, sc[0][0], sc[0][1], sc[0][2],
                            sc[1][0], sc[1][1], sc[1][2],
                            sc[2][0], sc[2][1], sc[2][2],
                            tc[0][0], tc[0][1],
                            tc[1][0], tc[1][1],
                            tc[2][0], tc[2][1],
                            intensity, zbuffer, model);
    }
    free(zbuffer);
    tgaFlipVertically(image);    
}

int main(int argc, char **argv){
    if (argc != 4){
        fprintf(stderr, "Usage: %s <objfile> <diffusemap> <outfile>\n", argv[0]);
        return -1;
    }

    Model *model = NULL;
    model = loadFromObj(argv[1]);
    
    if (!model){
        perror("loadFromObj");
        freeModel(model);
        return -1;
    }
    
    if (!loadDiffuseMap(model, argv[2])){
        perror("loadDiffuseMap");
        freeModel(model);
        return -1;
    }
    
    tgaFlipVertically(model->diffuse_map);
    tgaFlipHorizontally(model->diffuse_map);

    tgaImage *image = NULL;
    image = tgaNewImage(800, 800, RGB);

    if (!image){
        perror("tgaNewImage");
        freeModel(model);
        tgaFreeImage(image);
        return -1;
    }
    
    render(image, model);

    if (-1 == tgaSaveToFile(image, argv[3])){
        perror("tgaSateToFile");
        freeModel(model);
        tgaFreeImage(image);
        return -1;
    }

    freeModel(model);
    tgaFreeImage(image);

    return 0;
}
