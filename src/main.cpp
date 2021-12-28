#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <liquid.hpp>
#include <iostream>


#define WIDTH 400
#define HEIGHT 400

#define I(X,Y) (Y*WIDTH + X)*4

using namespace sf;

void draw_fluid(FluidCube* cube, Uint8* pixels);


int main()
{
    RenderWindow window(VideoMode(WIDTH, HEIGHT), "Yes Fluids");
    Vector2i ppos;

    Uint8 pixels[WIDTH*HEIGHT*4];
    FluidCube* my_cube = FluidCubeCreate(400, 1, 1, 1);
    
    for (int x=180; x<220; x++ )
    {
        for (int y=180; y<220; y++ )
        {
            FluidCubeAddDensity(my_cube, x, y, 100);
        }
    }

    while (window.isOpen())
    {
        Event event;
        while (window.pollEvent(event))
        {
            if (event.type == Event::Closed)
                window.close();
        }

        if(Mouse::isButtonPressed(Mouse::Left))
        {
            Vector2i pos = Mouse::getPosition(window);
            printf("x: %i, y: %i\n", pos.x, pos.y);
            if(pos.x > 0 && pos.x <= WIDTH && pos.y > 0 && pos.y <= HEIGHT )
            {
                FluidCubeAddDensity(my_cube, pos.x, pos.y, 100);
                FluidCubeAddVelocity(my_cube, pos.x, pos.y, pos.x - ppos.x, pos.y - ppos.y);
            }

            ppos = pos;
            
        }

        FluidCubeStep(my_cube);
        // now print it 
        draw_fluid(my_cube, pixels);
        Image img;
        img.create(WIDTH, HEIGHT, pixels);

        Texture texture;
        texture.create(WIDTH,HEIGHT);
        texture.loadFromImage(img);

        Sprite sprite;
        sprite.setTexture(texture);
        window.clear();
        window.draw(sprite);
        window.display();
    }

    return 0;
}

void draw_fluid(FluidCube* cube, Uint8* pixels)
{
    for(uint x = 0; x < (uint) cube->size; x++)
    {
        for(uint y = 0; y < (uint) cube->size; y++)
        {
            pixels[I(x,y)] = (Uint8) cube->density[y*cube->size + x]; // R
            pixels[I(x,y)+1] = (Uint8) cube->density[y*cube->size + x]; // G
            pixels[I(x,y)+2] = (Uint8) cube->density[y*cube->size + x]; // B
            pixels[I(x,y)+3] = (Uint8) 255; // A
        }
    }
}

void clamp(int* value, int min, int max)
{
    if(*value < min)
        *value = min;
    else if(*value > max)
        *value = max;
}