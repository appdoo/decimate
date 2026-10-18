/**
 * TargaImage.h
 *
 * Headerfile for the TargaImage class
 *
 * @author Andor Greissl
 * @creation: 27.09.2002
 */
#ifndef TARGAIMAGE_H
#define TARGAIMAGE_H

// C++ Includes
#ifdef WIN32
#include <windows.h>
#endif
#include <gl\gl.h>
#include <gl\glu.h>
#include <gl\glaux.h>
#include <fstream>

// Local Includes
#include "../engine/ext_gl.h"


/**
 * TGA FILE-Header
 *
 * Information taken from 
 * "http://www.organicbit.com/closecombat/formats/tga.html"
 * - thanks a lot....
 *
 * @author Andor Greissl
 * @creation: 27.09.2002
 */

#pragma pack(1)
typedef struct
{
    unsigned char  identsize;          // size of ID field that follows 18 byte header (0 usually)
    unsigned char  colormaptype;       // type of color map 0=none, 1=has palette
    unsigned char  imagetype;          // type of image 0=none,1=indexed,2=rgb,3=grey,+8=rle packed

    unsigned short colormapstart;      // first color map entry in palette
    unsigned short colormaplength;     // number of colors in palette
    unsigned char  colormapbits;       // number of bits per palette entry 15,16,24,32

    short          xstart;             // image x origin
    short          ystart;             // image y origin
    short          width;              // image width in pixels
    short          height;             // image height in pixels
    unsigned char  bits;               // image bits per pixel 8,16,24,32
    unsigned char  descriptor;         // image descriptor bits (vh flip bits)
    
    // pixel data follows header
    
} TGA_HEADER;
#pragma pack()



/**
 * Texture class
 *
 * @author Andor Greissl
 * @creation: 27.09.2002
 */

class TargaImage
{
private:
	// tga header
	TGA_HEADER tga_header;

	// texture width
	int width;
	// texture height
	int height;
	// bytes per pixel
	int bpp;
	// the image data
	void *image;

public:

	/**
	 * construct a new, empty Texture 
	 *
	 *
	 * @author Andor Greissl
	 * @creation: 07.11.2002
	 */
	TargaImage();

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
   TargaImage(int width, int height, int bpp, void *image);

	/**
	 * clean up Texture object
	 *
	 *
	 * @author Andor Greissl
	 * @creation: 07.11.2002
	 */
	~TargaImage();

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
   TargaImage *readTga(char *fileName);

	/**
	 * returns the rgb value for the location
	 * u, v
	 *
	 * @author Andor Greissl
	 * @creation: 07.11.2002
	 */
	void getNormal(float *n, float u, float v);


private:
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
   void correctTga(TGA_HEADER &tga_header, char *image);


   /**
    * reads an compressed TGA file
    *
    * @param in, reference to the file input stream
    * @param tga_header, reference to the tga_header
    *
    * @author Andor Greissl
    * @creation: 07.11.2002
    */
   TargaImage *readCompressedTga(std::ifstream &in, TGA_HEADER &tga_header);

   /**
   * reads an compressed TGA file
   *
   * @param in, reference to the file input stream
   * @param tga_header, reference to the tga_header
   *
   * @author Andor Greissl
   * @creation: 07.11.2002
   */
    TargaImage *readTga(std::ifstream &in, TGA_HEADER &tga_header);

};

#endif