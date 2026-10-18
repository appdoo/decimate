/***************************************************************************
 * Class to split a polygon into patches
 *
 * Functionality
 *
 * Splitter class can be constructed from
 * a vertex array and a face - indices array, to keep
 * the interface small and general
 *
 * After finding the patches, All Patch Elements can be
 * loaded by calling an iterator method
 *
 *
 * Description
 *
 * The splitter class builds a winged edge data structure
 * from the given face - indices list
 *
 * for each edge, a sharpness criteria is calculated and the 
 * edges are sorted according to their sharpness criteria
 **************************************************************************/
#ifndef SPLITTER_H
#define SPLITTER_H

// Local Includes
#include "..\math\math3D.h"
#include "..\math\heap.h"
#include "..\math\halfedge.h"
#include "chart.h"


class Splitter
{
   // number of halfedges
   Halfedge *halfedges;
   int      numHalfedges;
   // number of edges
   Edge  *edges;
   int   numEdges;
   // number of faces
   Face  *faces;
   int   numFaces;
   // number of points
   Vec3  *points;
   int   numPoints;

   // neighbour tags
   bool *neighbours;
   // visited tags
   bool *visited;
   // detected features
   bool *features;

public:
   /**
    * Constructs a Splitter class from the given
    * coordinates and face indices
    */
   Splitter	(Vec3 *points, int numPoints, Face *faces, int numFaces);

   /**
    * Splitter cleanup
    */
   ~Splitter();

private:
   /***************************************************************************
    * initializes the Splitter data structure and creates
    * a winged edge structure
    **************************************************************************/
   void store_as_winged(Vec3 *points, int numPoints, Face *faces, int numFaces);

   /***************************************************************************
    * calculate the sharpness for a halfedge
    **************************************************************************/
   float sharpness(Halfedge *w, Vec3 *normals);

   /***************************************************************************
   * calculate the sharpness for all halfedges
   **************************************************************************/
   float *sharpness();
   
   /***************************************************************************
    * calculate the sharpness for a string of halfedges
    **************************************************************************/
   float Splitter::sharpness(halfedge_array *s, float *criteria);

public:
  /***************************************************************************
   * returns a sharpness visualiuation consisting
   * out of edges and color, representing the sharpness criteria
   **************************************************************************/
   void sharpnessVisualization(int *num, Edge **edges, rgbaf **rgbvalues);
   
   /******************************************************************************
   * Markieren einer Ecke
   *****************************************************************************/
   void mark(Halfedge *w, bool mark);

   /******************************************************************************
   * Aussage darüber, ob eine Ecke markiert ist oder nicht
   *****************************************************************************/
   bool is_marked(Halfedge *w);

   /******************************************************************************
   * Aussage darüber, ob eine Ecke bereits als Nachbar markiert ist
   *****************************************************************************/
   bool is_neighbour(Halfedge *w);

   /******************************************************************************
    * Tiefensuche über die Winged - Edge - Struktur
    *****************************************************************************/
   halfedge_array *depth_first(Halfedge *start, float treshold, float *criteria, int maxlen, Vec3 *dir);

   /******************************************************************************
    * tags all elements of vector as features
    *****************************************************************************/
   void tag_features(halfedge_array *detected_features, bool *features);
   
   /******************************************************************************
    * tags all element-neighbours of vector as feature_neighbours
    *****************************************************************************/
   void tag_feature_neighbours(halfedge_array *detected_features);

   /******************************************************************************
   * expand a feature curve
   *****************************************************************************/
   void expand_feature_curve(
      Halfedge *start,
      int max_string_length,
      int min_feature_length,
      float treshold,
      float *criteria,
      bool *features
      );

   /***************************************************************************
    * calculate a treshold for the edges to use
    **************************************************************************/
   float calculateTreshold(float percentage, const float *criteria);

   /***************************************************************************
   * static method for the qsort functionality
   **************************************************************************/
   static int floatSorter(const void *a, const void *b);

   /***************************************************************************
   * detect the features
   **************************************************************************/
   bool *detect_features(float percentage, int min_feature_length, int max_string_length);

   /***************************************************************************
   * returns a sharpness visualiuation consisting
   * out of edges and color, representing the sharpness criteria
   **************************************************************************/
   void featureVisualization(
      float percentage, int min_feature_length, int max_string_length,
      int *num, Edge **edges);

   /***************************************************************************
   * the distance_to_features function
   * calculates a distance to the feature curves for each face
   **************************************************************************/
   float *distance_to_features(bool *features);
   float *Splitter::distance_to_features2(bool *features);

   /***************************************************************************
   * the distance_to_features_visualization function
   * returns an array of edges, containing only those edges
   * whose incident faces have a different distance_to_feature value
   **************************************************************************/
   void distance_to_features_visualization(
         float percentage, int min_feature_length, int max_string_length,
         int *num, Edge **edges
      );

   /***************************************************************************
   * returns a sharpness visualiuation consisting
   * out of edges and color, representing the sharpness criteria
   **************************************************************************/
   void seedVisualization(
      float percentage, int min_feature_length, int max_string_length,
      int *numEdges, Edge **edges, rgbf **rgbvalues);

   /******************************************************************************
    * return the distance for charts to be merged
    *****************************************************************************/
   float getMergeDistance(float *distances);

   /******************************************************************************
    * checks if distance_to_features of f is a local maximum
    *****************************************************************************/
   bool is_local_max(int face_index, float *distances);

   void remove_non_extremal(Halfedge *start, Chart *chart, halfedge_set *boundaries);
   bool is_linked(Halfedge *w, halfedge_set *boundaries);

   /******************************************************************************
    * merges the opposite chart (opp) with chart c
    * corrects the references and chart boundaries
    *****************************************************************************/
   void merge_chart(Chart *c, Chart *opp, Chart **references, halfedge_set *chart_boundaries);

   /******************************************************************************
    * chart expansion algorithm
    *****************************************************************************/
   std::set<Chart *> *expand_charts(float percentage, int min_feature_length, int max_string_length);

   /******************************************************************************
   * creates a chart visualization and splits the model along the points
   *****************************************************************************/
   void chart_visualization(
         float percentage, int min_feature_length, int max_string_length,
         int *numEdges, Edge **edges, rgbf **rgb_values
      );

   /******************************************************************************
   * sets a color value for each different chart
   *****************************************************************************/
   void colorize(std::set<Chart *> *charts, int *num, Edge **edges, rgbf **rgb_values);

   /******************************************************************************
    * Computes the seeds for a set of charts
    *****************************************************************************/
   int Splitter::compute_seed(Chart *chart);

   /******************************************************************************
   * Compute the seeds for a set of charts
   *****************************************************************************/
   int_set *compute_seeds(std::set<Chart *> *charts);

   /******************************************************************************
   * Calculates the cost of an edge
   * Sander et al.
   *****************************************************************************/
   double cost(Chart *chart, Face *f, Face *fopp);

   /******************************************************************************
   * Calculates the cost of an edge
   * Sander et al.
   *****************************************************************************/
   double cost(Face *f, Face *fopp);

   /******************************************************************************
    * check if this edge is linked with an incoming edge from
    * boundaries
    *****************************************************************************/
   bool is_linked_incoming(Halfedge * w, halfedge_set *boundaries);


   /******************************************************************************
    * check if this edge is linked with an outgoing edge from
    * boundaries
    *****************************************************************************/
   bool is_linked_outgoing(Halfedge * w, halfedge_set *boundaries);

   void remove_non_extremal2(Halfedge *start, Chart *chart, bool *boundaries);

   /******************************************************************************
    * expand a set of charts simultaneously
    * using an edge cost, containing face to chart - normal and 
    * planar distance
    *****************************************************************************/
   chart_set *expand_charts(chart_set *chart_alloc, int_set *seeds);

   /******************************************************************************
   * Different Approach for chart growing
   *
   *****************************************************************************/
   chart_set *lloyd_expand_charts();

   /******************************************************************************
   * initialize a boolean array
   *****************************************************************************/
   void initialize(bool *array, int size, bool value);

   float Splitter::distance(Halfedge *w);

   /******************************************************************************
    * calculates maximum distance
    *****************************************************************************/
   float mid_halfedge(Halfedge **mid_halfedge, Chart *chart);

   /******************************************************************************
   * calculates maximum distance
   *****************************************************************************/
   float max_dist(Halfedge **max, Halfedge *start, Chart *chart);

   /******************************************************************************
    * calculates chart diameter
    *****************************************************************************/
   float diameter(int *p1, int *p2, Chart *chart);

   /******************************************************************************
    * calculates chart geometric diameter
    * maximum geometric distance between two poits
    *****************************************************************************/
   float geodiameter(int *p1, int *p2, Chart *chart);

   /******************************************************************************
    * returns the surface area of the model
    *****************************************************************************/
   float surfacearea();

   /******************************************************************************
    * creates a chart visualization and splits the model along the points
    *****************************************************************************/
   void texture_visualization(
      float percentage, int min_feature_length, int max_string_length,
      int *numPoints, int *numFaces, Vec3 **points, Face **faces, Tex2 **texture, rgb **rgbvalues
   );
};


#endif
