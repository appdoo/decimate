/**
 * Texture.cpp
 *
 * Implementation of the Texture class
 *
 * the Texture class combines image loading and converting
 * with opengl - texture binding, mipmap filtering, bump mapping
 * etc...
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
#include <stdio.h>
#include "targaimage.h"



/**
 * Texture class implementation
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */

/**
 * construct a new, empty Texture 
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
TargaImage::TargaImage()
{
	TargaImage(0, 0, 0, NULL);
}


/**
 * construct a TargaImage object from the given parameters
 *
 * @param width, new texture width
 * @param height, new texture height
 * @param bpp, new texture bits per pixel
 * @param image, new image data
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
TargaImage::TargaImage(int width, int height, int bpp, void *image)
{
	// copy parameters
	this->width  = width;
	this->height = height;
	this->bpp    = bpp;
	this->image  = image;
}


/**
 * clean up TargaImage object
 *
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
TargaImage::~TargaImage()
{
	if (image) delete[] image;
}


/**
 * reads a tga image, and returns
 * a Texture object
 *
 * @param char*, name of the tga image
 *
 * TODO throw exceptions on fail
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
TargaImage *TargaImage::readTga(char *fileName)
{
   std::ifstream in(fileName, std::ios_base::binary);
	TGA_HEADER tga_header;

	if (!in)
		return NULL;

	// read tga header
	if (!in.read((char*)&tga_header, sizeof(TGA_HEADER)))
		return NULL;

	// can't handle color maps... (yet)
	if (tga_header.colormaptype != 0)
		return NULL;

	// process only 24/32 bpp (others not implemented yet)
	if (tga_header.bits!=24  &&  tga_header.bits!=32)
		return NULL;

	// skip image information field
	//if(!in.seekg(tga_header.identsize, ios::cur))
	//	return NULL;

	// check if tga is compressed, branch
	if (tga_header.imagetype & 0x08b)
	{
		return readTga(in, tga_header);
	}
	else
	{
		return readCompressedTga(in, tga_header);
	}
	return NULL;
}

/**
 * correct tga specific image data
 * there are faster ways to do this, but 
 * there is no need for speed so i left the code in 
 * a readable form
 *
 * @param tga_header, a reference to the TGA Header of the image
 * @param image, a char pointer to the image data
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
void TargaImage::correctTga(TGA_HEADER &tga_header, char *image)
{
	// calculate bytes per pixel
	int byteSize = tga_header.bits >> 3;

	// calculate used memory
	int memsize = tga_header.width * tga_header.height * byteSize;

	// get start and end positions
	unsigned char *lp  = (unsigned char *) image;
	unsigned char *end = lp + memsize;


	// declare swap variables
	unsigned char r, g, b;

	while (lp < end)
	{
		// get tga specified bgr
		b = lp[0];
		g = lp[1];
		r = lp[2];

		// write std rgb
		lp[0] = r;
		lp[1] = g;
		lp[2] = b;

		// increment pointer
		lp += byteSize;
	}
}


/**
 * reads an compressed TGA file
 *
 * @param in, reference to the file input stream
 * @param tga_header, reference to the tga_header
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
TargaImage *TargaImage::readTga(std::ifstream &in, TGA_HEADER &tga_header)
{
	// allocate space for image data
	int memsize = tga_header.width * tga_header.height * (tga_header.bits >> 3);
	char *image = new char[memsize];

	// load image data
	if (!in.read(image, memsize))
	{
		delete[] image;
		return NULL;
	}
	
	// correct tga specific image data
	correctTga(tga_header, image);

	// return a new texture object
	return new TargaImage(tga_header.width, tga_header.height, tga_header.bits, image);
}


/**
 * reads an compressed TGA file
 *
 * @param in, reference to the file input stream
 * @param tga_header, reference to the tga_header
 *
 * @author Andor Greissl
 * @creation: 07.11.2002
 */
TargaImage *TargaImage::readCompressedTga(std::ifstream &in, TGA_HEADER &tga_header)
{
	return NULL;
}
