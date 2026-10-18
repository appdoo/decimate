// ============================================================================
//
// Definition der Schnittstelle zum Normal - Sampling
//
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 4 $
// 
// ============================================================================
#ifndef TEXTURESAMPLE_H
#define TEXTURESAMPLE_H

// Local Includes
#include "../math/math3D.h"
#include "../math/vec.h"

namespace normalmap {


/******************************************************************************
 * computes normal values for texture map
 *****************************************************************************/
void mapnormals(
                   int texSize,
                   rgb *mem,
                   int  numLowPoins,
                   Vec3 *lowPoints,
                   Vec3 *lowNormals,
                   Tex2 *lowTex,
                   int   numLowFaces, 
                   Face *lowFaces,

                   Vec3 *highPoints,
                   Vec3 *highNormals,
                   int  numHighFaces,
                   Face *highFaces
                   );

} // end namespace

#endif
