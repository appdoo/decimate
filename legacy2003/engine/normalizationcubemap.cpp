/**
 * NormalizationCubeMap.cpp
 *
 *
 * Implementation of the NormalizationCubeMap class
 *
 * NormalizationCubeMap can be used to normalize vectors within
 * the RegisterCombiners
 *
 */
#include <math.h>
#include "normalizationcubemap.h"




/******************************************************************************
 * constructor for a normalization cube map
 *****************************************************************************/
NormalizationCubeMap::NormalizationCubeMap(int size)
{
	unsigned char *texture = new unsigned char[size * size * 3];
	float nrm[3];

	glGenTextures(1, &texturename);
	glBindTexture(GL_TEXTURE_CUBE_MAP_EXT, texturename);

  glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP_EXT, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  for (int i = 0; i < 6; i++) {
    for (int y = 0; y < size; y++) {
      for (int x = 0; x < size; x++) {
        getNormal(nrm, i, x, y, size);
				//
				rangeCompress(nrm);
				//
        texture[(y * size + x) * 3 + 0] = (unsigned char) nrm[0];
				texture[(y * size + x) * 3 + 1] = (unsigned char) nrm[1];
        texture[(y * size + x) * 3 + 2] = (unsigned char) nrm[2];
      }
    }
    glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT + i, 0,	GL_RGB8, size, size, 0, GL_RGB, GL_UNSIGNED_BYTE, texture);
	}
	delete[] texture;
}


/******************************************************************************
 * destructor
 *****************************************************************************/
NormalizationCubeMap::~NormalizationCubeMap()
{
	glDeleteTextures(1, &texturename);
}


/******************************************************************************
 * do a range compress
 *****************************************************************************/
void NormalizationCubeMap::rangeCompress(float *n)
{
	n[0] = (n[0] + 1.0f) * 127.0f;
	n[1] = (n[1] + 1.0f) * 127.0f;
	n[2] = (n[2] + 1.0f) * 127.0f;
}


/******************************************************************************
 * bring a vector to unit length
 *****************************************************************************/
void NormalizationCubeMap::norm(float *nrm)
{
	// normalize the notmal vector
  float len = (float) (1.0 / sqrt(nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2]));

	// bring vector to length 1
  nrm[0] *= len;
  nrm[1] *= len;
  nrm[2] *= len;
}


/******************************************************************************
 * create a normal the cube map
 *****************************************************************************/
void NormalizationCubeMap::getNormal(float *n, int faceIndex, int x, int y, int size)
{
  float w = (float)x / (float) (size - 1);
  float h = (float)y / (float) (size - 1);
  w = w + w - 1.0f;
  h = h + h - 1.0f;

  switch (faceIndex) {
		case 0:
			n[0] = 1.0f;
			n[1] = -h;
			n[2] = -w;
			break;
		case 1:
			n[0] = -1.0f;
			n[1] = -h;
			n[2] = w;
			break;
		case 2:
			n[0] = w;
			n[1] = 1.0f;
			n[2] = h;
			break;
		case 3:
			n[0] = w;
			n[1] = -1.0f;
			n[2] = -h;
			break;
		case 4:
			n[0] = w;
			n[1] = -h;
			n[2] = 1.0f;
			break;
		case 5:
			n[0] = -w;
			n[1] = -h;
			n[2] = -1.0f;
			break;
		default:
			break;
  }
	norm(n);
}


/******************************************************************************
 * open gl bind texture
 *****************************************************************************/
void NormalizationCubeMap::bind()
{
	glBindTexture(GL_TEXTURE_CUBE_MAP_EXT, texturename);
}
