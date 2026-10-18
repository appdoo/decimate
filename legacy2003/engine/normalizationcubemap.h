/**
 * NormalizationCubeMap.h
 *
 *
 * HeaderFile for the NormalizationCubeMap class
 *
 *
 * NormalizationCubeMap can be used to normalize vectors within
 * the RegisterCombiners
 *
 */

#ifndef NORMALIZATIONCUBEMAP_H
#define NORMALIZATIONCUBEMAP_H

#include "ext_gl.h"


class NormalizationCubeMap
{
public:
	GLuint texturename;


public:
	/******************************************************************************
	 * constructor for a normalization cube map
	 *****************************************************************************/
	NormalizationCubeMap(int size);


	/******************************************************************************
	 * destructor
	 *****************************************************************************/
	~NormalizationCubeMap();


	/******************************************************************************
	 * open gl bind texture
	 *****************************************************************************/
	void bind();

private:
	/******************************************************************************
	 * bring a vector to unit length
	 *****************************************************************************/
	 void norm(float *nrm);


	/******************************************************************************
	 * do a range compress
	 *****************************************************************************/
	void rangeCompress(float *nrm);


	/******************************************************************************
	 * create a normal for the cube map
	 *****************************************************************************/
	void getNormal(float *nrm, int faceIndex, int x, int y, int size);
};




#endif