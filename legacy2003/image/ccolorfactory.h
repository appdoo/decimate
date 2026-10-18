//=============================================================================
//
// CColorFactory : Class Definition
//
// Stellt Farben für die Visualisierung bereit
//
// $Author: Andor $
// $Date: 6.12.03 1:14 $
// $Revision: 1 $
// 
//=============================================================================
#ifndef CCOLORFACORY_H
#define CCOLORFACORY_H

class CColorFactory
{
public:
   // initialize factory
   static void Initialize()
   {
      srand(6870967);
   }

   // creates a color
   static rgbf HtmlColor()
   {
      // temporary color
      rgbf col;

      // generate a html compatible color
      col.r = (rand() * 4.0f * 0x33 + 0x33 ) / ((float)RAND_MAX * 255.0f);
      col.g = (rand() * 4.0f * 0x33 + 0x33 ) / ((float)RAND_MAX * 255.0f);
      col.b = (rand() * 4.0f * 0x33 + 0x33 ) / ((float)RAND_MAX * 255.0f);

      return col;
   }
};




#endif