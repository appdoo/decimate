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
 *
 * $Author: Andor $
 * $Date: 17.12.03 2:05 $
 * $Revision: 7 $
 *
 ***************************************************************************/

// Local Includes
#include "splitter.h"
#include "chart_packing.h"


/***************************************************************************
 * Constructs a Splitter class from the given
 * coordinates and face indices
 **************************************************************************/
Splitter::Splitter(Vec3 *points, int numPoints, Face *faces, int numFaces)
{
   // initialize pointers
   this->halfedges = NULL;
   this->edges  = NULL;
   this->faces  = NULL;
   this->points = NULL;
   this->numHalfedges = 0;
   this->numEdges  = 0;
   this->numFaces  = 0;
   this->numPoints = 0;
   // setup data structure
   store_as_winged(points, numPoints, faces, numFaces);
}

/***************************************************************************
 * Splitter cleanup
 **************************************************************************/
Splitter::~Splitter()
{
   if (halfedges) delete[] halfedges;
   if (edges)  delete[] edges;
   if (faces)  delete[] faces;
   if (points) delete[] points;
}

/***************************************************************************
 * initializes the Splitter data structure and creates
 * a winged edge structure
 **************************************************************************/
void Splitter::store_as_winged(Vec3 *points, int numPoints, Face *faces, int numFaces)
{
   // copy points
   if (points && numPoints>0) {
      this->points = new Vec3[numPoints];
      this->numPoints = numPoints;
      //
      memcpy(this->points, points, sizeof(Vec3) * numPoints);
   }
   // copy faces
   if (faces && numFaces > 0) {
      this->faces = new Face[numFaces];
      this->numFaces = numFaces;
      //
      memcpy(this->faces, faces, sizeof(Face) * numFaces);
   }
   // create a winged edge structure
   HalfEdgeMesh::Create(numPoints, numFaces, faces, &numHalfedges, &halfedges);
}

/***************************************************************************
 * calculate the sharpness for a halfedge
 **************************************************************************/
float Splitter::sharpness(Halfedge *w, Vec3 *normals)
{
   // check if halfedge is a border
   // and has no oppositehalfedge
   // in this case the maximum sharpness of 1.0 is returned
   if (!w->sym)
      return 1.0f;
   // otherwise, calculate angle between face normals
   // as a sharpness criteria
   return 1.0f - fabs(normals[w->face] * normals[w->sym->face]);
}

/***************************************************************************
 * calculate the sharpness for all halfedges
 **************************************************************************/
float *Splitter::sharpness()
{
   // allocate space for face normals
   Vec3 *normals = new Vec3[numFaces];
   float *criteria = new float[numHalfedges];
   // get all face normals
   for (int i=0; i<numFaces; i++) {
      normals[i] = normal(&faces[i], points);
   }
   // get sharpness criteria for all halfedges
   for (int i=0; i<numHalfedges; i++) {
      criteria[i] = sharpness(&halfedges[i], normals);
   }
   // cleanup and exit
   delete[] normals;
   return criteria;
}

/***************************************************************************
 * calculate the sharpness for a string of halfedges
 **************************************************************************/
float Splitter::sharpness(halfedge_array *s, float *criteria)
{
   // validate vector
   if (!s) return 0.0f;
   // calculate sum
   float sum = 0.0f;
   for (halfedge_array::iterator i=s->begin(); i!=s->end(); i++)
      sum += criteria[(*i) - halfedges];
   return sum;
}

/***************************************************************************
 * returns a sharpness visualiuation consisting
 * out of edges and color, representing the sharpness criteria
 **************************************************************************/
void Splitter::sharpnessVisualization(int *num, Edge **edges, rgbaf **rgbvalues)
{
   // calculate sharpness for all halfedges
   float *criteria = sharpness();
   // create a marker for visited edges
   bool *marker = new bool[numHalfedges];
   memset(marker, false, sizeof(bool) * numHalfedges);
   Edge *tmp_edge = new Edge[numHalfedges];
   *rgbvalues = new rgbaf[numPoints];
   for (int i=0; i<numPoints; i++) {
      (*rgbvalues)[i].r = 0.9f;
      (*rgbvalues)[i].g = 0.9f;
      (*rgbvalues)[i].b = 0.9f;
      (*rgbvalues)[i].a = 0.5f;
   }
   int numEdges = 0;
   // convert from winged edge to edge
   for (int i=0; i<numHalfedges; i++) {
      // check marker
      if (marker[i] == true)
         continue;
      // create a edge from winged edge
      int a = halfedges[i].end;
      int b = halfedges[i].next->next->end;
      tmp_edge[numEdges].a = a;
      tmp_edge[numEdges].b = b;
      float rgba = MIN(1.0f - criteria[i], (*rgbvalues)[a].r);
      (*rgbvalues)[a].r = rgba;
      (*rgbvalues)[a].g = rgba;
      (*rgbvalues)[a].a = 1.0f - rgba;
      //(*rgbvalues)[a].b = rgba;
      //(*rgbvalues)[a].a = MAX(criteria[i], (*rgbvalues)[a].a);
      float rgbb = MIN(1.0f - criteria[i], (*rgbvalues)[b].r);
      (*rgbvalues)[b].r = rgbb;
      (*rgbvalues)[b].g = rgbb;
      (*rgbvalues)[b].a = 1.0f - rgbb;
      //(*rgbvalues)[b].b = rgbb;
      //(*rgbvalues)[b].a = MAX(criteria[i], (*rgbvalues)[b].a);
      numEdges++;
      marker[i]=true;
      // set symmetric marker
      if (halfedges[i].sym) {
         int symidx = halfedges[i].sym - &halfedges[0];
         marker[symidx]=true;
      }
   }
   // reallocate
   *edges = new Edge[numEdges];
   *num = numEdges;
   // copy results
   memcpy(*edges, tmp_edge, sizeof(Edge) * numEdges);
   // delete temporary arrays
   delete[] tmp_edge;
   delete[] criteria;
}

/******************************************************************************
 * Markieren einer Ecke
 *****************************************************************************/
void Splitter::mark(Halfedge *w, bool mark)
{
   if (w) visited[w->end] = mark;
}

/******************************************************************************
 * Aussage darüber, ob eine Ecke markiert ist oder nicht
 *****************************************************************************/
bool Splitter::is_marked(Halfedge *w)
{
   return w && visited[w->end];
}

/******************************************************************************
 * Aussage darüber, ob eine Ecke bereits als Nachbar markiert ist
 *****************************************************************************/
bool Splitter::is_neighbour(Halfedge *w)
{
   return w && neighbours[w - halfedges] || (w && w->sym && neighbours[w->sym - halfedges]);
}

/******************************************************************************
 * Tiefensuche über die Winged - Edge - Struktur
 *****************************************************************************/
halfedge_array *Splitter::depth_first(Halfedge *start, float treshold, float *criteria, int maxlen, Vec3 *dir)
{
   // define a result, and a maximum vector
   halfedge_array *max = NULL;
   float max_sharpness = 0.0f;

   // mark wegdes as currently visited
   mark(start, true);

   // get the first outgoing halfedge
   Halfedge *loop = start->next;
   while (maxlen>1 && loop && loop->sym!=start && start->end==loop->next->next->end) {
      // calculate search direction
      Vec3 search_vec = (points[loop->end] - points[start->end]).norm();

      // no halfedge of path is allowed to be a neighbour
      if (!((search_vec * (*dir))<0.0f || is_neighbour(loop) || is_marked(loop))) {
         // find the part with maximum sharpness
         halfedge_array *result = depth_first(loop, treshold, criteria, maxlen - 1, dir);
         // check for the first result
         if (max == NULL) {
            max = result;
            max_sharpness = sharpness(max, criteria);
         }
         // take the one with higher sharpness
         else if (sharpness(result, criteria) > max_sharpness) {
            delete max;
            max = result;
            max_sharpness = sharpness(max, criteria);
         }
         else {
            // cleanup result allocation
            delete result;
         }
      }
      // go to the next outgoing halfedge (if exists)
      if (loop->sym) loop = loop->sym->next;
      else           loop = NULL;
   }
   // delete marks
   mark(start, false);
   // check if we have a maximum
   if (!max) {
      max = new halfedge_array;
      max->push_back(start);
   }
   // insert first halfedge
   else max->insert(max->begin(), start);
   // returns the result
   return max;
}

/******************************************************************************
 * tags all elements of vector as features
 *****************************************************************************/
void Splitter::tag_features(halfedge_array *s, bool *features)
{
   for (halfedge_array::iterator i=s->begin(); i!=s->end(); i++) {
      features[(*i) - halfedges] = true;
   }
}

/******************************************************************************
 * tags all element-neighbours of vector as feature_neighbours
 *****************************************************************************/
void Splitter::tag_feature_neighbours(std::vector<Halfedge *> *s)
{
   for (halfedge_array::iterator i=s->begin(); i!=s->end(); i++) {
      // get halfedge
      Halfedge *start = *i;
      Halfedge *loop  = start;
      // loop around the endpoint
      do {
         // marke halfedge
         neighbours[loop - halfedges] = true;
         neighbours[loop->next - halfedges] = true;
         neighbours[loop->next->next - halfedges] = true;
         // marke symmetric halfedge
         if (loop->sym) {
            neighbours[loop->sym - halfedges] = true;
            neighbours[loop->sym->next - halfedges] = true;
            neighbours[loop->sym->next->next - halfedges] = true;
         }
         // go to next edge around endpoint
         loop = loop->next->sym;
      } while (loop && loop!=start);
   }
}

/******************************************************************************
 * expand a feature curve
 *****************************************************************************/
void Splitter::expand_feature_curve(
   Halfedge *start,
   int max_string_length,
   int min_feature_length,
   float treshold,
   float *criteria,
   bool *features
   )
{
   //our detected feature curve
   halfedge_array detected_feature;
   // loop over halfedge and symmetric halfedge
   for (int i=0; i<2; i++) {
      // check start condition
      if (!start) continue;
      // treat symmetric halfedge, halfedge
      Halfedge *h = start;
      start = start->sym;
      // detect feature
      detected_feature.clear();
      //detected_feature.push_back(h);
      float break_condition = 0.0f;
      do {
         break_condition = 0.0f;
         // calculate starting direction vector
         Vec3 direction = points[h->end] - points[h->next->next->end];
         direction.norm();
         // search next string, sarting with edge h
         halfedge_array *s = depth_first(h, treshold, criteria, max_string_length, &direction);
         int string_length = s->size();
         if (s && s->size() > 1) {
            mark(h, true);
            h = (*s)[1];
            detected_feature.push_back(h);
            break_condition = sharpness(s, criteria);
         }
         if (s) delete s;
      }
      while (break_condition > max_string_length * treshold);
      // remove the temporary marks
      //for (halfedge_array::iterator i=detected_feature.begin(); i!=detected_feature.end(); i++)
         //mark(*i, false);
      // check for min_feature_length
      if ((int)detected_feature.size() > min_feature_length) {
         tag_features(&detected_feature, features);
         tag_feature_neighbours(&detected_feature);
      }
   }
}

/***************************************************************************
 * calculate a treshold for the edges to use
 **************************************************************************/
float Splitter::calculateTreshold(float percentage, const float *criteria)
{
   float *dummy = new float[numHalfedges];
   memcpy(dummy, criteria, sizeof(float) * numHalfedges);
   qsort(dummy, numHalfedges, sizeof(float), Splitter::floatSorter);
   //
   // check borders
   if (percentage < 0.01f) percentage = 0.01f;
   if (percentage > 1.0f) percentage = 1.0f;
   // get index
   int index = (int) ((numHalfedges - 1) * percentage);
   float treshold = dummy[index];
   delete[] dummy;
   //
   //return treshold value
   return treshold;
}

/***************************************************************************
 * static method for the qsort functionality
 **************************************************************************/
int Splitter::floatSorter(const void *a, const void *b)
{
   // sort descending
   return SGN(*(float*)b - *(float*)a);
}

/***************************************************************************
 * detect the features
 **************************************************************************/
bool *Splitter::detect_features(float percentage, int min_feature_length, int max_string_length)
{
   // calculate sharpness for all halfedges
   float *criteria = sharpness();
   // allocate markers
   neighbours = new bool[numHalfedges];
   visited = new bool[numPoints];
   features = new bool[numHalfedges];
   //
   // initialize markers
   for (int i=0; i<numHalfedges; i++) {
      // set marker to false
      features[i]   = false;
      neighbours[i] = false;
   }
   // mark points for depth - first search
   for (int i=0; i<numPoints; i++) {
      // set marker to false
      visited[i]    = false;
   }

   // recalculate the treshold from a percentage of edges to use
   // sharpness criteria
   float treshold = calculateTreshold(percentage, criteria);
   
   int num = 0;

   // go through all halfedges finding sharpness contours
   printf("Searching for features...\n");
   for (int i=0; i<numHalfedges; i++) {
      // check treshold
      if (criteria[i] < treshold)
         continue;
      // check neighbourhood
      if (is_neighbour(&halfedges[i]) || is_neighbour(halfedges[i].sym))
         continue;
      // check feature
      if (features[i])
         continue;
      // use halfedge as starting point
      expand_feature_curve(&halfedges[i], max_string_length, min_feature_length, treshold, criteria, features);
      //features[i] = true;
      num++;
   }
   printf("Used halfedges [%d] of [%d]\n", num, numHalfedges);

   // cleanup with the sharpness criteria
   delete[] criteria;
   delete[] visited;
   delete[] neighbours;
   // return result
   return features;
}

/***************************************************************************
 * returns a feature visualiuation consisting
 * out of edges and color, representing the features
 **************************************************************************/
void Splitter::featureVisualization (
      float percentage, int min_feature_length, int max_string_length,
      int *num, Edge **edges
   )
{
   // calculate sharpness for all halfedges
   bool *features = detect_features(percentage, min_feature_length, max_string_length);
   // create a marker for visited edges
   bool *marker = new bool[numHalfedges];
   memset(marker, false, sizeof(bool) * numHalfedges);
   Edge *tmp_edge = new Edge[numHalfedges];

   int numEdges = 0;
   // convert from winged edge to edge
   for (int i=0; i<numHalfedges; i++) {
      // check marker
      if (marker[i] == true || !features[i])
         continue;
      // create a edge from winged edge
      tmp_edge[numEdges].a = halfedges[i].end;
      tmp_edge[numEdges].b = halfedges[i].next->next->end;
      numEdges++;
      marker[i]=true;
      // set symmetric marker
      if (halfedges[i].sym)
         marker[halfedges[i].sym - halfedges] = true;
   }
   // reallocate
   *edges = new Edge[numEdges];
   *num = numEdges;
   // copy results
   memcpy(*edges, tmp_edge, sizeof(Edge) * numEdges);
   // delete temporary arrays
   delete[] tmp_edge;
   delete[] features;
}

/***************************************************************************
 * the distance_to_features function
 * calculates a distance to the feature curves for each face
 **************************************************************************/
float *Splitter::distance_to_features(bool *features)
{
   // allocate space for distance information
   float *distances = new float[numFaces];
   // initialize distance to negative value
   for (int i=0; i<numFaces; i++) {
      distances[i] = -1.0f;
   }
   // put all the feature curves into an heap
   Heap<Halfedge *, float> front;
   for (int i=0; i<numHalfedges; i++) {
      // check if halfedge is a feature
      if (features[i]) {
         front.push(&halfedges[i], 0.0f);
         if (halfedges[i].sym) front.push(halfedges[i].sym, 0.0f);
      }
   }

   // propagate the front
   Halfedge *start;
   float value;
   // process all heap - elements
   while (front.pop(&start, &value)) {
      // Step 1: enter face distances
      Halfedge *loop  = start;
      // loop around the endpoint
      do {
         // check if face has the correct distance
         if (distances[loop->face]==-1 || distances[loop->face] > value) {
            // set face distance
            distances[loop->face] = value;
         }
         // go to next edge around endpoint
         if (loop->next && loop->next->sym) {
            loop = loop->next->sym;
         }
         else loop = NULL;
      } while (loop && loop!=start);
      // increment value
      value = value + 0.01f;

      // Step 2: write back edges
      loop  = start;
      // loop around the endpoint
      do {
         // loop inside the triangle
         Halfedge *w = loop;
         for (int i=0; i<3; i++) {
            // next side of the face
            w = w->next;
            if (w->sym && (distances[w->sym->face]==-1 || distances[w->sym->face]>value)) {
               // add halfedge to heap
               front.push(w->sym, value);
            }
         }
         // go to next edge around endpoint
         if (loop->next && loop->next->sym) {
            loop = loop->next->sym;
         }
         else loop = NULL;
      } while (loop && loop!=start);
   }
   // return the information
   return distances;
}

/***************************************************************************
 * the distance_to_features function
 * calculates a distance to the feature curves for each face
 **************************************************************************/
float *Splitter::distance_to_features2(bool *features)
{
   // allocate space for distance information
   float *distances = new float[numFaces];
   // initialize distance to negative value
   for (int i=0; i<numFaces; i++) {
      distances[i] = -1.0f;
   }
   // put all the feature curves into an heap
   Heap<Halfedge *, float> front;
   for (int i=0; i<numHalfedges; i++) {
      // check if halfedge is a feature
      if (features[i] || !halfedges[i].sym) {
         front.push(&halfedges[i], 0.0f);
         if (halfedges[i].sym) front.push(halfedges[i].sym, 0.0f);
      }
   }

   // propagate the front
   Halfedge *w;
   float value;
   // process all heap - elements
   while (front.pop(&w, &value)) {
      // check if there is anything todo for this edge
      if (distances[w->face]!=-1 /*&& distances[w->face] <= value*/) {
         distances[w->face]=MIN(value, distances[w->face]);
         continue;
      }
      // enter the face distance_to_feature information
      distances[w->face] = value;
      // calculate some distance information
      Vec3 a = (points[faces[w->face].b] + points[faces[w->face].a] + points[faces[w->face].c]) * 0.333f;
      a -= points[faces[w->face].a];
      // increment the distance
      value = value + 0.01f;//a.length(); //0.01f;
      // check opposite triangles 1
      for (int i=0; i<2; i++) {
         w = w->next;
         if (w->sym && (distances[w->sym->face]==-1.0f || distances[w->sym->face]>value)) {
            // add halfedge to heap
            front.push(w->sym, value);
         }
      }
   }
   // return the information
   return distances;
}

/***************************************************************************
 * returns a sharpness visualiuation consisting
 * out of edges and color, representing the sharpness criteria
 **************************************************************************/
void Splitter::seedVisualization(
   float percentage, int min_feature_length, int max_string_length,
   int *numEdges, Edge **edges, rgbf **rgbvalues)
{
   featureVisualization(percentage, min_feature_length, max_string_length, numEdges, edges);
   // calculate sharpness for all halfedges
   bool *features = detect_features(percentage, min_feature_length, max_string_length);
   // get distance_to_feature values
   printf("calculating distance to features\n");
   float *distances = distance_to_features2(features);
   // allocate colors
   *rgbvalues = new rgbf[numPoints];
   //initialize to white
   for (int i=0; i<numPoints; i++) {
      (*rgbvalues)[i].r = 1.0f;
      (*rgbvalues)[i].g = 1.0f;
      (*rgbvalues)[i].b = 1.0f;
   }
   for (int i=0; i<numFaces; i++) {
      // check for local maximum
      if (is_local_max(i, distances)) {
         (*rgbvalues)[faces[i].a].r = 1.0f;
         (*rgbvalues)[faces[i].a].g = 0.0f;
         (*rgbvalues)[faces[i].a].b = 0.0f;
         (*rgbvalues)[faces[i].b].r = 1.0f;
         (*rgbvalues)[faces[i].b].g = 0.0f;
         (*rgbvalues)[faces[i].b].b = 0.0f;
         (*rgbvalues)[faces[i].c].r = 1.0f;
         (*rgbvalues)[faces[i].c].g = 0.0f;
         (*rgbvalues)[faces[i].c].b = 0.0f;      }
   }
   delete[] features;
   delete[] distances;
}

/***************************************************************************
 * the distance_to_features_visualization function
 * returns an array of edges, containing only those edges
 * whose incident faces have a different distance_to_feature value
 **************************************************************************/
void Splitter::distance_to_features_visualization(
      float percentage, int min_feature_length, int max_string_length,
      int *num, Edge **edges
   )
{
   // calculate sharpness for all halfedges
   bool *features = detect_features(percentage, min_feature_length, max_string_length);
   // get distance_to_feature values
   printf("calculating distance to features\n");
   float *distances = distance_to_features(features);
   printf("done...\n");
   // create a marker for visited edges
   bool *marker = new bool[numHalfedges];
   memset(marker, false, sizeof(bool) * numHalfedges);
   Edge *tmp_edge = new Edge[numHalfedges];

   int numEdges = 0;
   // convert from winged edge to edge
   for (int i=0; i<numHalfedges; i++) {
      // check marker
      if (marker[i] == true)
         continue;
      // check if we have an edge with two faces
      if (!halfedges[i].sym)
         continue;
      // check equi-distance criteria
      if (distances[halfedges[i].face] == distances[halfedges[i].sym->face])
         continue;
      // create a edge from winged edge
      tmp_edge[numEdges].a = halfedges[i].end;
      tmp_edge[numEdges].b = halfedges[i].next->next->end;
      numEdges++;
      // set marker
      marker[i]=true;
      // set symmetric marker
      if (halfedges[i].sym) {
         marker[halfedges[i].sym - halfedges] = true;
      }
   }
   // reallocate
   *edges = new Edge[numEdges];
   *num = numEdges;
   // copy results
   memcpy(*edges, tmp_edge, sizeof(Edge) * numEdges);
   // cleanup 
   delete[] tmp_edge;
   delete[] features;
   delete[] distances;
}

/******************************************************************************
 * return the distance for charts to be merged
 *****************************************************************************/
float Splitter::getMergeDistance(float *distances)
{
   // start with first value for maximum
   float max = distances[0];
   // search for max within the rest
   for (int i=1; i<numFaces; i++) {
      max = MAX(max, distances[i]);
   }
   // return an certain amount of distance
   return max / 4.0f;
}

/******************************************************************************
 * checks if distance_to_features of f is a local maximum
 *
 * we say a face has maximum distance, when all incident
 * faces of all face - points have a lower or eqal distance
 *****************************************************************************/
bool Splitter::is_local_max(int face_index, float *distances)
{
   // get face distance
   float d = distances[face_index];
   if (d<=0.0f)
      return false;
   // get corresponding halfedge
   Halfedge *w = &halfedges[face_index * 3];
   // check distances of corresponding faces
   for (int i=0; i<3; i++) {
      // check if distance of f is greater or equal
      Halfedge *loop = w->next->sym;
      while (loop && loop!=w) {
         if (d < distances[loop->face])
            return false;
         for (int j=0; j<3; j++) {
            if (loop->sym && loop->sym->face!=face_index && d < distances[loop->sym->face])
               return false;
            loop = loop->next;
         }
         // go to next edge around endpoint
         loop = loop->next->sym;
      }
      // go to next halfedge
      w = w->next;
   }
   // we have a local maximum
   return true;
}


bool Splitter::is_linked(Halfedge *w, halfedge_set *boundaries)
{
   // check for valid halfedge
   if (!w) 
      return true;
   // go to next outgoing halfedge
   Halfedge *loop = w->next;
   while (loop && loop->sym!=w) {
      // check if there is any halfedge linked
      if (boundaries->find(loop) != boundaries->end())
         return true;
      //if (loop->sym && boundaries->find(loop->sym) != boundaries->end())
      //   return true;
      // go to next edge around endpoint
      if (loop->sym) loop = loop->sym->next;
      else           loop = NULL;
   }
   return false;
}

void Splitter::remove_non_extremal(Halfedge *start, Chart *chart, halfedge_set *boundaries)
{
   for (int i=0; i<3;i++) {
      start = start->next;
      if (!(is_linked(start, boundaries) && is_linked(start->sym, boundaries))) {
         boundaries->erase(start);
         boundaries->erase(start->sym);
         chart->erase(start, points);
         chart->erase(start->sym, points);
      }
   }
/*
   halfedge *loop = start;
   for (int i=0; i<3; i++) {
      // check if chart_boundaries contains halfedge
      if (loop->sym && chart->contains(loop->sym->face)) {
         boundaries->erase(loop);
         boundaries->erase(loop->sym);
         chart->erase(start);
         chart->erase(start->sym);
      }
      loop = loop->next;
   }
   */
}
/******************************************************************************
 * merges the opposite chart (opp) with chart c
 * corrects the references and chart boundaries
 *****************************************************************************/
void Splitter::merge_chart(Chart *c, Chart *opp, Chart **references, halfedge_set *chart_boundaries)
{
   // for all faces belonging to opposite chart
   for (Chart::face_iterator it=opp->begin_face(); it!=opp->end_face(); ++it) {
      // get the halfedge index
      Halfedge *w = &halfedges[(*it) * 3];
      // for all three face edges do
      for (int j=0; j<3; j++) {
         // erase chart boundaries belonging to chart c and opp
         if (w->sym && c->contains(w->sym->face)) {
            chart_boundaries->erase(w);
            chart_boundaries->erase(w->sym);
         }
         w = w->next;
      }
      // set reference from opp to c
      references[*it] = c;
   }
   // merge chart data
   c->merge(opp);
}

/******************************************************************************
 * chart expansion algorithm
 *****************************************************************************/
std::set<Chart *> *Splitter::expand_charts(float percentage, int min_feature_length, int max_string_length)
{
   // detect features
   bool *features = detect_features(percentage, min_feature_length, max_string_length);
   // get distance to features for each face
   float *distances = distance_to_features2(features);
   delete[] features;
   float epsilon = getMergeDistance(distances);
   // Heap sorted by dist(facet(halfedge))
   Heap<Halfedge *, float> heap;
   // face to chart index chart
   Chart **charts = new Chart*[numFaces];
   // charts
   std::set<Chart *> *chart_alloc = new std::set<Chart *>;
   // Initialize Heap
   for (int i=0; i<numFaces; i++) {
      // clear chart pointer
      charts[i] = NULL;
      // check if distance_to_features is a local maximum
      if (is_local_max(i, distances)) {
         // allocate a new chart
         charts[i] = new Chart(i, distances[i]);
         chart_alloc->insert(charts[i]);
         // add the halfefdges of F to Heap
         Halfedge *w = &halfedges[i * 3];
         heap.push(w,             -distances[i]);
         heap.push(w->next,       -distances[i]);
         heap.push(w->next->next, -distances[i]);
      }
   }
   // initialized with all the edges of the surface
   std::set<Halfedge *> chart_boundaries;
   for (int i=0; i<numHalfedges; i++) {
      chart_boundaries.insert(&halfedges[i]);
   }
   // local variables for chart growing
   Halfedge *h = NULL;
   float value = 0.0f;
   // Charts-growing phase
   while(heap.pop(&h, &value)) {
      // check if we have a symmetric halfedge
      if (h->sym == NULL)
         continue;
      // continue with chart growing
      int f = h->face;
      int fopp = h->sym->face;
      // check if we have a chart for fopp
      if (charts[fopp] == NULL) {
         // add fopp to chart(F)
         charts[f]->add(fopp, distances[fopp]);
         charts[fopp] = charts[f];
         // remove E from chart_boundaries
         // for halfedges we have to erase both halfedge and symmetric halfedge
         chart_boundaries.erase(h);
         chart_boundaries.erase(h->sym);
         // remove non-extremal edges from chart_boundaries,
         // (i.e. edges that do not link two other chart boundary edges)
         remove_non_extremal(h->sym, charts[f], &chart_boundaries);

         // add the halfedges of Fopp belonging to
         // chart_boundaries to Heap
         Halfedge *loop = h->sym;
         for (int i=0; i<3; i++) {
            // check if chart_boundaries contains halfedge
            if (chart_boundaries.find(loop)!=chart_boundaries.end()) {
               heap.push(loop, -distances[loop->face]);
            }
            loop = loop->next;
         }
      }
      // check if to charts should be merged
      else if ((charts[fopp] != charts[f]) &&
               (charts[f]->max_dist() - distances[f] < epsilon) &&
               (charts[fopp]->max_dist() - distances[f] < epsilon )) 
      {
         //merge chart(F) and chart(Fopp)
         Chart *copp = charts[fopp];
         merge_chart(charts[f], charts[fopp], charts, &chart_boundaries);
         chart_alloc->erase(copp);
         delete copp;
         // print debug information
         printf("%d\n", chart_alloc->size());
      }
   }
   // enter boundary information into chart
   for (halfedge_set::iterator it=chart_boundaries.begin(); it!=chart_boundaries.end(); ++it) {
      if (charts[(*it)->face])
         charts[(*it)->face]->add(*it, points);
      else printf("lost face %d\n", (*it)->face);
   }
   // cleanup temporary stuff
   delete[] distances;
   delete[] charts;
   // return a chart array
   return chart_alloc;
}

/******************************************************************************
 * creates a chart visualization and splits the model along the points
 *****************************************************************************/
void Splitter::chart_visualization(
      float percentage, int min_feature_length, int max_string_length,
      int *num, Edge **edges, rgbf **rgb_values
   )
{
   // get the charts
   //chart_set *charts = expand_charts(percentage, min_feature_length, max_string_length);
   chart_set *charts = lloyd_expand_charts();
   // initialize color
   colorize(charts, num, edges, rgb_values);
   // delete chart set
   for (std::set<Chart *>::iterator it=charts->begin(); it!=charts->end(); it++) {
      delete *it;
   }
   delete charts;
}


/******************************************************************************
 * sets a color value for each different chart
 *****************************************************************************/
void Splitter::colorize(std::set<Chart *> *charts, int *num, Edge **edges, rgbf **rgb_values)
{
   // initialize color
   *rgb_values = new rgbf[numPoints];
   memset(*rgb_values, 0, sizeof(rgbf) * numPoints);
   rgbf col;
   srand(6870967);

   // create a marker for visited edges
   bool *marker = new bool[numHalfedges];
   memset(marker, false, sizeof(bool) * numHalfedges);
   Edge *tmp_edge = new Edge[numHalfedges];
   int numEdges = 0;

   for (std::set<Chart *>::iterator it=charts->begin(); it!=charts->end(); it++) {
      // process chart
      Chart *chart = *it;
      // generate a html compatible color
      col.r = (rand() * 4.0f * 0x33 + 0x33 ) / ((float)RAND_MAX * 255.0f);
      col.g = (rand() * 4.0f * 0x33 + 0x33 ) / ((float)RAND_MAX * 255.0f);
      col.b = (rand() * 4.0f * 0x33 + 0x33 ) / ((float)RAND_MAX * 255.0f);
      // process all faces
      for (Chart::face_iterator fit=chart->begin_face(); fit!=chart->end_face(); ++fit) {
         // load face
         Face *f = &faces[*fit];
         // set color
         (*rgb_values)[f->a] = col;
         (*rgb_values)[f->b] = col;
         (*rgb_values)[f->c] = col;
      }
      // process all halfedges
      for (Chart::halfedge_iterator it = chart->begin_halfedge(); it!=chart->end_halfedge(); it++) {
         // get halfedge
         Halfedge *w = *it;
         int i = w - halfedges;
         // check marker
         if (marker[i] == true)
            continue;
         // create a edge from winged edge
         tmp_edge[numEdges].a = halfedges[i].end;
         tmp_edge[numEdges].b = halfedges[i].next->next->end;
         numEdges++;
         // set marker
         marker[i]=true;
         // set symmetric marker
         if (halfedges[i].sym) {
            marker[halfedges[i].sym - halfedges] = true;
         }
      }
   }

   // reallocate
   *edges = new Edge[numEdges];
   *num = numEdges;
   // copy results
   memcpy(*edges, tmp_edge, sizeof(Edge) * numEdges);
}

/******************************************************************************
 * Compute the seeds for a set of charts
 *****************************************************************************/
int Splitter::compute_seed(Chart *chart)
{
   Halfedge *mid;
   mid_halfedge(&mid, chart);
   return mid->face;
}

/******************************************************************************
 * Compute the seeds for a set of charts
 *****************************************************************************/
std::set<int> *Splitter::compute_seeds(std::set<Chart *> *charts)
{
   // allocate seed set
   std::set<int> *seeds = new std::set<int>;
   // iterate over all charts
   for (std::set<Chart *>::iterator it=charts->begin(); it!=charts->end(); it++) {
      // compute seed for the face
      int seed = compute_seed(*it);
      // at seed to set of seeds
      seeds->insert(seed);
   }
   return seeds;
}

/******************************************************************************
 * Calculates the cost of an edge
 * Sander et al.
 *****************************************************************************/
double Splitter::cost(Chart *chart, Face *f, Face *fopp)
{
   // get normals and centroids
   Vec3 fnrm    = chart->average_normal();
   Vec3 fmid    = (points[f->a] + points[f->b] + points[f->c]) / 3.0f;
   Vec3 foppnrm = normal(fopp, points);
   Vec3 foppmid = (points[fopp->a] + points[fopp->b] + points[fopp->c]) / 3.0f;
   //
   float area = chart->area();
   float perimeter = chart->perimeter();
   area += normal(fopp, points).length() * 0.5f;
   Halfedge *w = &halfedges[(fopp - faces) * 3];
   for (int i=0; i<3; i++) {
      w = w->next;
      float length = distance(w);
      if (chart->chart_boundaries.find(w->sym) == chart->chart_boundaries.end())
         perimeter += length;
      else
         perimeter -= length;
   }
   float ratio = (perimeter * perimeter) / (4.0f * (float)M_PI * area);
   return (1.0f - (fnrm * foppnrm)) * (fmid - foppmid).length() + ratio;
}

/******************************************************************************
 * Calculates the cost of an edge
 * Sander et al.
 *****************************************************************************/
double Splitter::cost(Face *f, Face *fopp)
{
   // get normals and centroids
   Vec3 fnrm    = normal(f, points);
   Vec3 fmid    = (points[f->a] + points[f->b] + points[f->c]) / 3.0f;
   Vec3 foppnrm = normal(fopp, points);
   Vec3 foppmid = (points[fopp->a] + points[fopp->b] + points[fopp->c]) / 3.0f;
   //
   return (1.0 - (fnrm * foppnrm)) * (fmid - foppmid).length();
}


/******************************************************************************
 * check if this edge is linked with an incoming edge from
 * boundaries
 *****************************************************************************/
bool Splitter::is_linked_incoming(Halfedge * w, halfedge_set *boundaries)
{
   Halfedge *loop = w;
   do {
      // go to incoming halfedge
      loop = loop->next->next;
      // check if halfedge belongs to boundaries
      if (loop!=w->sym && loop!=w->next->next && boundaries->find(loop) != boundaries->end())
         return true;
      // go to face
      loop = loop->sym;
   } while (loop && loop!=w);
   return false;
}

/******************************************************************************
 * check if this edge is linked with an outgoing edge from
 * boundaries
 *****************************************************************************/
bool Splitter::is_linked_outgoing(Halfedge * w, halfedge_set *boundaries)
{
   Halfedge *loop = w;
   do {
      // go to incoming halfedge
      loop = loop->next;
      // check if halfedge belongs to boundaries
      if (loop!=w->sym && boundaries->find(loop) != boundaries->end())
         return true;
      // go to face
      loop = loop->sym;
   } while (loop && loop!=w);
   return false;
}

void Splitter::remove_non_extremal2(Halfedge *start, Chart *chart, bool *boundaries)
{
   for (int i=0; i<3;i++) {
      start = start->next;
      if (start->sym && chart->contains(start->sym->face)) {
         boundaries[start - halfedges] = false;
         boundaries[start->sym - halfedges] = false;
         chart->erase(start->sym, points);
         //chart->erase(start, points);
      }
   }
}

/******************************************************************************
 * expand a set of charts simultaneously
 * using an edge cost, containing face to chart - normal and 
 * planar distance
 *****************************************************************************/
chart_set * Splitter::expand_charts(chart_set *chart_alloc, int_set *seeds)
{
   // reset chart allocation
   chart_alloc->clear();
   // Heap sorted by cost(halfedge)
   Heap<Halfedge *, double> heap;
   // face to chart index chart
   Chart **charts = new Chart*[numFaces];
   memset(charts, 0, sizeof(Chart *) * numFaces);
   // Initialize Heap
   for (int_set::iterator it=seeds->begin(); it!=seeds->end(); it++) {
      // get seed information and create chart
      int seed_face = *it;
      // create a Chart for each seed
      Chart *chart = new Chart(seed_face, normal(&faces[seed_face], points));
      charts[seed_face] = chart;
      chart_alloc->insert(chart);

      // get the halfedges belonging to the face
      Halfedge *w = &halfedges[seed_face * 3];
      // push the three halfedges onto heap
      for (int i=0; i<3; i++) {
         // add to chart boundaries
         chart->add(w, points);
         // add halfedge to heap for search
         // if we have an opposite face
         if (w->sym) {
            heap.push(w, cost(&faces[*it], &faces[w->sym->face]));
         }
         // loop inside face
         w = w->next;
      }
   }
   seeds->clear();
   // initialized with all the edges of the surface
   bool *chart_boundaries = new bool[numHalfedges];
   for (int i=0; i<numHalfedges; i++) {
      chart_boundaries[i] = true;
   }
   // local variables for chart growing
   Halfedge *h = NULL;
   double value = 0.0f;
   // Charts-growing phase
   while(heap.pop(&h, &value)) {
      // check if we have a symmetric halfedge
      if (h->sym == NULL)
         continue;
      // continue with chart growing
      int f    = h->face;
      int fopp = h->sym->face;
      // check if we have a chart for fopp
      if (charts[fopp] == NULL) {
         // add fopp to chart(F)
         charts[f]->add(fopp, normal(&faces[fopp], points));
         charts[fopp] = charts[f];
         // remove E from chart_boundaries
         // for halfedges we have to erase both halfedge and symmetric halfedge
         chart_boundaries[h - halfedges] = false;
         chart_boundaries[h->sym - halfedges] = false;
         charts[f]->erase(h, points);
         charts[f]->erase(h->sym, points);
         // remove non-extremal edges from chart_boundaries,
         // (i.e. edges that do not link two other chart boundary edges)
         remove_non_extremal2(h->sym, charts[f], chart_boundaries);

         // add the halfedges of Fopp belonging to
         // chart_boundaries to Heap
         Halfedge *loop = h->sym->next;
         for (int i=0; i<2; i++) {
            // check if chart_boundaries contains halfedge
            if (chart_boundaries[loop - halfedges]) {
               charts[f]->add(loop, points);
               if (loop->sym) {
                  double cst = cost(charts[f], &faces[loop->face], &faces[loop->sym->face]);
                  heap.push(loop, cst);
               }
            }

            loop = loop->next;
         }
      }
   }
/*   // create boundary structure
   for (int i=0; i<numHalfedges; i++) {
      if (chart_boundaries[i] && charts[halfedges[i].face]) {
         charts[halfedges[i].face]->add(&halfedges[i]);
      }
   }*/
   // cleanup temporary stuff
   delete[] charts;
   delete[] chart_boundaries;
   // return a chart array
   return chart_alloc;
}

/*
chart_set * Splitter::expand_charts(chart_set *chart_alloc, int_set *seeds)
{
   // reset chart allocation
   chart_alloc->clear();
   // Heap sorted by cost(halfedge)
   Heap<Halfedge *, float> heap;
   // face to chart index chart
   Chart **charts = new Chart*[numFaces];
   memset(charts, 0, sizeof(Chart *) * numFaces);
   // Initialize Heap
   for (int_set::iterator it=seeds->begin(); it!=seeds->end(); it++) {
      // create a Chart for each seed
      charts[*it] = new Chart(*it, normal(&faces[*it]));
      chart_alloc->insert(charts[*it]);

      // get the halfedges belonging to the face
      Halfedge *w = &halfedges[*it * 3];
      // push the three halfedges onto heap
      for (int i=0; i<3; i++) {
         // loop inside face
         w = w->next;
         // add to chart borders
         charts[*it]->add(w);
      }
      push_minimal_cost(&heap, *it, charts);
   }
   // initialized with all the edges of the surface
   halfedge_set chart_boundaries;
   for (int i=0; i<numHalfedges; i++) {
      chart_boundaries.insert(&halfedges[i]);
   }
   // local variables for chart growing
   Halfedge *h = NULL;
   float value = 0.0f;
   // Charts-growing phase
   while(heap.pop(&h, &value)) {
      // check if we have a symmetric halfedge
      if (h->sym == NULL)
         continue;
      // continue with chart growing
      int f = h->face;
      int fopp = h->sym->face;
      // check if we have a chart for fopp
      if (charts[fopp] == NULL) {
         // add fopp to chart(F)
         charts[f]->add(fopp, normal(&faces[fopp]));
         charts[fopp] = charts[f];
         // remove E from chart_boundaries
         // for halfedges we have to erase both halfedge and symmetric halfedge
         chart_boundaries.erase(h);
         chart_boundaries.erase(h->sym);
         charts[f]->erase(h);
         charts[f]->erase(h->sym);
         // remove non-extremal edges from chart_boundaries,
         // (i.e. edges that do not link two other chart boundary edges)
         remove_non_extremal2(h->sym, charts[f], &chart_boundaries);

         // add the halfedges of Fopp belonging to
         // chart_boundaries to Heap
         Halfedge *loop = h->sym;
         for (int i=0; i<3; i++) {
            // check if chart_boundaries contains halfedge
            if (chart_boundaries.find(loop)!=chart_boundaries.end()) {
               //if (loop->sym)
               //   heap.push(loop, cost(charts[f], &faces[loop->face], &faces[loop->sym->face]));
               charts[f]->add(loop);
            }
            loop = loop->next;
         }
         push_minimal_cost(&heap, f, charts);
      }
   }
   // cleanup temporary stuff
   delete[] charts;
   // return a chart array
   return chart_alloc;
}
*/

/******************************************************************************
 * initialize a boolean array
 *****************************************************************************/
void Splitter::initialize(bool *array, int size, bool value)
{
   for (int i=0; i<size; i++)
      array[i] = value;
}

/******************************************************************************
 * calculates maximum distance
 *****************************************************************************/
float Splitter::mid_halfedge(Halfedge **mid_halfedge, Chart *chart)
{
   // heap for algorithm
   Heap<Halfedge *, float> heap;
   visited = new bool[numPoints];
   initialize(visited, numPoints, false);
   Halfedge *start, *loop;
   float max   = 0.0f;
   float value = 0.0f;

   for (halfedge_iterator it=chart->begin_halfedge(); it!=chart->end_halfedge(); it++) {
      start = *it;
      // mark starting point
      mark(start, true);
      heap.push(start, 1.0f);
      // initialize heap
      *mid_halfedge  = start;

      loop = start->next;
      while (loop && loop->sym!=start) {
         // check if chart contains halfedge
         if (!is_marked(loop) && (chart->contains(loop->face)))
            heap.push(loop, 1.0f);
         // next outgoing
         if (loop->sym) loop = loop->sym->next;
         else           loop = NULL;
      }
   }
   // apply search
   while (heap.pop(&start, &value)) {
      // check for new maximum
      if (value > max) {
         max = value;
         *mid_halfedge = start;
      }
      // already been there
      if (is_marked(start))
         continue;
      // mark point as visited
      mark(start, true);
      loop = start->next;
      while (loop && loop->sym!=start) {
         // check if chart contains halfedge
         if (!is_marked(loop) && (chart->contains(loop->face)))
            heap.push(loop, value + 1.0f);
         // next outgoing
         if (loop->sym) loop = loop->sym->next;
         else           loop = NULL;
      }
   }
   delete[] visited;
   return max;
}

/******************************************************************************
 * returns the length of a halfedge
 *****************************************************************************/
float Splitter::distance(Halfedge *w)
{
   return (points[w->end] - points[w->next->next->end]).length();
}

/******************************************************************************
 * calculates maximum distance
 *****************************************************************************/
float Splitter::max_dist(Halfedge **max_halfedge, Halfedge *start, Chart *chart)
{
   // heap for algorithm
   Heap<Halfedge *, float> heap; 
   initialize(visited, numPoints, false);

   // mark starting point
   mark(start->next->next, true);
   heap.push(start, 1.0f);
   // initialize heap
   Halfedge *loop = start->next->next->sym;
   *max_halfedge  = start;
   float max   = 1.0f;
   float value = 0.0f;
   while (loop && loop!=start) {
      // check if chart contains halfedge
      if (chart->contains(loop->face))
         heap.push(loop, distance(loop));
      // next outgoing
      loop = loop->next->next->sym;
   }

   // apply search
   while (heap.pop(&start, &value)) {
      // check for new maximum
      if (value > max) {
         max = value;
         *max_halfedge = start;
      }
      // already been there
      if (is_marked(start))
         continue;
      // mark point as visited
      mark(start, true);
      loop = start->next;
      while (loop && loop->sym!=start) {
         // check if chart contains halfedge
         if (!is_marked(loop) && (chart->contains(loop->face)))
            heap.push(loop, value + distance(loop));
         // next outgoing
         if (loop->sym) loop = loop->sym->next;
         else           loop = NULL;
      }
   }
   return max;
}

/******************************************************************************
 * calculates chart diameter
 *****************************************************************************/
float Splitter::diameter(int *p1, int *p2, Chart *chart)
{
   // define locals
   Halfedge *max_halfedge;
   float max = 0.0f;
   visited = new bool[numPoints];
   // calculate maximum distance for each chart
   for (halfedge_iterator it=chart->begin_halfedge(); it!=chart->end_halfedge(); it++) {
      // calculate maximum distance for each halfedge
      float dist = max_dist(&max_halfedge, *it, chart);
      // check for new maximum
      if (dist > max) {
         // set new maximum
         *p1 = (*it)->next->next->end;
         *p2 = max_halfedge->end;
         max = dist;
      }
   }
   // cleanup
   delete[] visited;
   return max;
}

/******************************************************************************
 * calculates chart geometric diameter
 * maximum geometric distance between two poits
 *****************************************************************************/
float Splitter::geodiameter(int *p1, int *p2, Chart *chart)
{
   // define locals
   Halfedge *max_halfedge;
   float max = 0.0f;
   visited = new bool[numPoints];
   // calculate maximum distance for each chart
   for (halfedge_iterator it=chart->begin_halfedge(); it!=chart->end_halfedge(); it++) {
      // calculate maximum distance for each halfedge
      for (halfedge_iterator jt=it; jt!=chart->end_halfedge(); jt++) {
         float dist = (points[(*it)->end] - points[(*jt)->end]).length();
         // check for new maximum
         if (dist > max) {
            // set new maximum
            *p1 = (*it)->end;
            *p2 = (*jt)->end;
            max = dist;
         }
      }
   }
   // cleanup
   delete[] visited;
   return max;
}

/******************************************************************************
 * Different Approach for chart growing
 *
 *****************************************************************************/
chart_set *Splitter::lloyd_expand_charts()
{
   // allocate Chart Array
   chart_set *charts = new chart_set;
   // seeds, start with face 0 (or a random face index)
   int_set *seeds = new int_set;
   int_set *new_seeds = NULL;
   for (int i=0; i<15; i++)
      seeds->insert(i);

   int iter=0;
   int maxiter = 15;
   // main texture loop
   while (++iter < maxiter)
   {
      // print seeds
      for (int_iterator it=seeds->begin(); it!=seeds->end(); it++)
         printf("%d ", *it);
      printf("\n");
      // cleanup
      for (chart_set::iterator it=charts->begin(); it!=charts->end(); it++)
         delete *it;
      charts = expand_charts(charts, seeds);
      new_seeds = compute_seeds(charts);
      delete seeds;
      seeds = new_seeds;
   }
   return charts;
}

/******************************************************************************
 * returns the surface area of the model
 *****************************************************************************/
float Splitter::surfacearea()
{
   float modelsize = 0.0f;

   for (int i=0; i<numFaces; i++) {
      modelsize += normal_raw(points, &faces[i]).length();
   }
   return 0.5f * modelsize;
}

/******************************************************************************
 * creates a chart visualization and splits the model along the points
 *****************************************************************************/
void Splitter::texture_visualization(
      float percentage, int min_feature_length, int max_string_length,
      int *numPoints, int *numFaces, Vec3 **points, Face **faces, Tex2 **texture,
      rgb **rgbvalues
   )
{
   // get the charts
   //chart_set *charts = expand_charts(percentage, min_feature_length, max_string_length);
   chart_set *charts = lloyd_expand_charts();

   int *path;
   int numPath;

   float surface = surfacearea();

   std::vector<Vec3 *> points_coll;
   std::vector<Tex2 *> tex_coll;
   std::vector<Face *> face_coll;

   std::vector<int> points_size;
   std::vector<int> face_size;

   std::vector<float> diameter;

   chartpacking::Horizon horizon(256);
   horizon.scale = 0.65f / sqrtf(surface) * horizon.resolution;

   *rgbvalues = new rgb[horizon.resolution * horizon.resolution];

   for (int i=0; i<(horizon.resolution * horizon.resolution); i++) {
      (*rgbvalues)[i].r = 128;
      (*rgbvalues)[i].g = 128;
      (*rgbvalues)[i].b = 255;
   }

   Vec3 *dstPoints   = NULL;
   int  numDstPoints = 0;
   Tex2 *dstTexture  = NULL;
   int  numDstFaces  = 0;
   Face *dstFaces    = NULL;

   // delete chart set
   int chart_index = 1;
   for (std::set<Chart *>::iterator it=charts->begin(); it!=charts->end(); it++) {
      Chart *chart = *it;
      int p1=0;
      int p2=1;

      printf("\n\nParametrizing chart\n");
      printf("diameter ");
      float dia = geodiameter(&p1, &p2, chart);
      printf("= %lf, %d, %d\n", dia, p1, p2);

      chart->parametrize_lscm(this->points, this->faces, p1, p2, numPoints, numFaces, points, faces, texture, &numPath, &path);
      chartpacking::chartpack(horizon, numPath / 2, (Edge *)path, *numPoints, *texture);

      // enter path
      for (int i=0; i<horizon.resolution; i++) {
         int offset = horizon.resolution * horizon.horizon[i] + i;
         if (offset<0 || offset>= horizon.resolution * horizon.resolution)
            continue;
         (*rgbvalues)[offset].r = 0xff;
         (*rgbvalues)[offset].g = 0x00;
         (*rgbvalues)[offset].b = 0x00;
      }

      Vec3 *tmpPoints   = new Vec3[numDstPoints + *numPoints];
      Tex2 *tmpTexture  = new Tex2[numDstPoints + *numPoints];
      Face *tmpFaces    = new Face[numDstFaces + *numFaces];

      // copy old stuff
      if (dstPoints) {
         memcpy(tmpPoints, dstPoints, sizeof(Vec3) * numDstPoints);
         memcpy(tmpTexture, dstTexture, sizeof(Tex2) * numDstPoints);
         memcpy(tmpFaces, dstFaces, sizeof(Face) * numDstFaces);

         delete[] dstPoints;
         delete[] dstTexture;
         delete[] dstFaces;
      }
      dstPoints  = tmpPoints;
      dstTexture = tmpTexture;
      dstFaces   = tmpFaces;

      // copy new stuff
      memcpy(&dstPoints[numDstPoints], *points, sizeof(Vec3) * *numPoints);
      memcpy(&dstTexture[numDstPoints], *texture, sizeof(Tex2) * *numPoints);

      for (int i=0; i<(*numFaces); i++) {
         dstFaces[numDstFaces + i].a = (*faces)[i].a + numDstPoints;
         dstFaces[numDstFaces + i].b = (*faces)[i].b + numDstPoints;
         dstFaces[numDstFaces + i].c = (*faces)[i].c + numDstPoints;
      }
      numDstFaces  += *numFaces;
      numDstPoints += *numPoints;

      // cleanup
      delete[] *points;
      delete[] *faces;
      delete[] path;
      delete[] *texture;

      delete *it;
   }
   delete charts;

   *numFaces = numDstFaces;
   *numPoints = numDstPoints;
   *points = dstPoints;
   *texture = dstTexture;
   *faces = dstFaces;
}


