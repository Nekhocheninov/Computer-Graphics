# Computer Graphics in C

Step-by-step development of a tiny software rasterizer based on [this material](https://github.com/ssloy/tinyrenderer/wiki/Lesson-0:-getting-started) in C.

### How To Start

1. Clone the repo:
    ```
    git clone https://github.com/Nekhocheninov/ComputerGraphics.git
    ```
2. Use [GNU make](https://gnuwin32.sourceforge.net/packages/make.htm) utility for building project:
    ```
    make -f makefile
    ```
3. Run the command to render:
    ```
    .\render obj\cat.obj obj\cat_diffuse.tga result.tga
    ```
# Sample

[Step 1](https://github.com/Nekhocheninov/ComputerGraphics/tree/wireframe-rendering): Creating a wireframe renderer.

<img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/wireframe-rendering/img_1.png" width="800">

[Step 2](https://github.com/Nekhocheninov/ComputerGraphics/tree/triangle): Filling triangles.

<img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/triangle/img_1.png" width="400"> <img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/triangle/img_2.png" width="400">

[Step 3](https://github.com/Nekhocheninov/ComputerGraphics/tree/z-buffer): Hidden faces removal.

<img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/z-buffer/img_1.png" width="800">

[Step 4](https://github.com/Nekhocheninov/ComputerGraphics/tree/texture): Adding texture.

<img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/texture/img_1.png" width="800">

[Step 5](https://github.com/Nekhocheninov/ComputerGraphics/tree/perspective-projection): Adding perspective projection.

<img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/perspective-projection/img_1.png" width="800">

[Step 6](https://github.com/Nekhocheninov/ComputerGraphics/tree/мoving-the-camera): Moving the camera.

Look! It has a tail!

<img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/мoving-the-camera/img_2.png" width="400"> <img src="https://github.com/Nekhocheninov/ComputerGraphics/blob/мoving-the-camera/img_1.png" width="400">
