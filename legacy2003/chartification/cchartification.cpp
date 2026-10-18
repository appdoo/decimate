// ============================================================================
//
// Chartification : Class Implementation
//
// Aufteilung in TexturStücke mittels Chart-Growing
//
// $Author: Andor $
// $Date: 17.12.03 2:05 $
// $Revision: 27 $
// 
// ============================================================================

// local includes
#include "../math/algorithm.h"
#include "../image/ccolorfactory.h"
#include "../utility.h"
#include "chartification.h"
#include "cchartboundary.h"
#include "cconformalmapper.h"
#include "catlaspacker.h"

#include "cchartset.h"

namespace chartification {
// =============================================================================
//
// CChartification Class Implementaion
//
//
// ----------------------------------------------------------------------------
// creates an CChartification object
// ----------------------------------------------------------------------------
CChartification::CChartification(int numPoints, Vec3 *points, int numFaces, Face *faces)
{
   this->numPoints = numPoints;
   this->numFaces  = numFaces;
   this->points    = new Vec3[numPoints];
   this->faces     = new Face[numFaces];

   // copy data
   memcpy(this->points, points, sizeof(Vec3) * numPoints);
   memcpy(this->faces, faces, sizeof(Face) * numFaces);

   // create a halfedge representation
   HalfEdgeMesh::Create(numPoints, numFaces, this->faces, &numHalfedges, &halfedges);

   // pre - calculate face normals and midpoints
   faceNormals     = new Vec3[numFaces];
   faceMidpoints   = new Vec3[numFaces];

   for (int i=0; i<numFaces; i++) {
      faceNormals[i]   = normal_raw(points, &faces[i]);
      faceMidpoints[i] = middle(points, &faces[i]);
   }
}


// ----------------------------------------------------------------------------
// deletes an CChartification object
// ----------------------------------------------------------------------------
CChartification::~CChartification()
{
   if (points)        delete[] points;
   if (faces)         delete[] faces;
   if (halfedges)     delete[] halfedges;
   if (faceNormals)   delete[] faceNormals;
   if (faceMidpoints) delete[] faceMidpoints;
}


// ----------------------------------------------------------------------------
// starts the chartification process
// ----------------------------------------------------------------------------
void CChartification::Chartify(int numDesiredCharts, CChartSet *charts)
{
   // create a chart representation
   chartMap      = new CChart*[numFaces];
   // Create Boundary structure
   boundary = new CChartBoundary(numPoints, numHalfedges, halfedges);

   // bootstrap process
   CChart *chart =  charts->Create(0);

   setup_timer();
   start_timer();

   // main chartification loop
   do {
      show_progress();
      ExpandCharts(charts);

      if (charts->Size() < numDesiredCharts) 
      {
         chart = charts->Create(lastface);
      }
   }
   while (!ComputeSeeds(charts, boundary));

   // check chartification for charts containing only one face
   // or charts whose irregulartity is to high
   while (charts->RemoveInvalidCharts())
   {
      ExpandCharts(charts);
   }
   while (Straighten(boundary));


   // create chart boundary structure
   CreateVertexStructure(charts, boundary);

   printf("\rAbzahl Charts %d, Millisekunden %lf\n\n", charts->Size(), stop_timer());
   delete[] chartMap;
   delete boundary;
}


// ----------------------------------------------------------------------------
// returns the cost for moving from one face to the other
// ----------------------------------------------------------------------------
double CChartification::EdgeCost(Halfedge *h)
{
   // get normals and centroids
   Vec3 nrm1 = faceNormals[h->face];
   Vec3 nrm2 = faceNormals[h->sym->face];

   // switch to chart normal
   if (chartMap[h->face]) nrm1 = chartMap[h->face]->Normal();
   else if (h->sym && chartMap[h->sym->face]) nrm2 = chartMap[h->sym->face]->Normal();

   // get Factor for difference in normals
   double nrmFactor = (double)(nrm1 * nrm2) / sqrt(nrm1.length2() * nrm2.length2());
   double distance  = Distance(h);

   return (1.0 - nrmFactor) * distance;
}


// ----------------------------------------------------------------------------
// distance between two faces
// ----------------------------------------------------------------------------
double CChartification::Distance(Halfedge *h)
{
   return (faceMidpoints[h->face] - faceMidpoints[h->sym->face]).length();
}


// ----------------------------------------------------------------------------
// Initialize heap and visited flags
// ----------------------------------------------------------------------------
void CChartification::InitializeSeedSearch(Heap<Halfedge *, double> *heap, bool *visited, CChartBoundary *boundary)
{
   // clear marker for visited faces
   for (int i=0; i<numFaces; i++) visited[i] = false;

   // initialize heap witj all boundary edges
   for (int i=0; i<numHalfedges; i++) {
      // process halfedge
      Halfedge *h = &halfedges[i];
      // check for boundary
      if (boundary->IsBoundary(h) && chartMap[h->face]) {
         // process boundary edge
         visited[h->face] = true;
         // enter seed position
         chartMap[h->face]->searchIndex = h->face;

         // push to neighbour halfedges on heap
         for (int i=0; i<2; i++) {
            // process the two neighbours
            h = h->next;
            // check if we have a valid ongoing halfedge of the chart
            if (h->sym && !visited[h->sym->face] && (chartMap[h->face] == chartMap[h->sym->face]))
               heap->push(h, Distance(h));
         }
      }
   }
}


// ----------------------------------------------------------------------------
// applies temorary seeds to seeds and checks for seed cycles
// ----------------------------------------------------------------------------
bool CChartification::ApplySeeds(CChartSet *charts)
{
   // flag to determine if theres a cycle within the seeds
   bool cycleDetected = true;

   // iterate over all charts
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // compute a new seed
      int seed = (*it)->searchIndex;

      // keep track of cycles in the seed
      if(seed>=0 && seed<numFaces) {
         cycleDetected &= (*it)->SetSeed(seed, faceNormals[seed]);
      }
      else cycleDetected = false;
   }
   return cycleDetected;
}


// ----------------------------------------------------------------------------
// Compute the seeds for a set of charts
// ----------------------------------------------------------------------------
bool CChartification::ComputeSeeds(CChartSet *charts, CChartBoundary *boundary)
{
   // heap for search algorithm
   Heap<Halfedge *, double> heap;
   Halfedge                 *loop;
   double                   value;

   // marker for already visited faces
   bool *visited = new bool[numFaces];

   // Initialize heap and visited flags
   InitializeSeedSearch(&heap, visited, boundary);

   // search the graph of the mesh
   while (heap.pop(&loop, &value)) {
      // switch to opposite halfedge
      loop = loop->sym;
      // test if weve been already there
      if (visited[loop->face]) continue;

      // mark face
      visited[loop->face] = true;

      // enter new maximum
      chartMap[loop->face]->searchIndex = loop->face;

      // push to neighbour halfedges on heap
      for (int i=0; i<2; i++) {
         // process the two neighbours
         loop = loop->next;
         // check if we have a valid ongoing halfedge of the chart
         if (loop->sym && !visited[loop->sym->face] && (chartMap[loop->face] == chartMap[loop->sym->face]))
            heap.push(loop, value + Distance(loop));
      }
   }
   delete[] visited;

   // reset seeds for minimal charts
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      if ((*it)->ClusterCount() < 4) (*it)->searchIndex = (*it)->Seed();
   }
   return ApplySeeds(charts);
}


// ----------------------------------------------------------------------------
// repair heap
// ----------------------------------------------------------------------------
void CChartification::RepairHeap(CChart *chart, Heap<Halfedge *, double> *heap)
{
   Halfedge *h;
   Halfedge *last = NULL;
   double    value;
   const int size = heap->size();

   // die momentane heap Implemntierung
   // setzt nur den Zähler auf 0, die Daten selbst werden
   // nicht verändert
   heap->clear();

   for (int i=0; i<size; i++) {
      // get value from original heap
      heap->get(i, &h, &value);

      if (chartMap[h->face] == chart) {
         // correct value
         value = EdgeCost(h);
      }
      heap->push(h, value);
   }
   //heap->repair();
}


// ----------------------------------------------------------------------------
// expand a set of charts simultaneously
// using an edge cost, containing face to chart - normal and 
// planar distance
// ----------------------------------------------------------------------------
void CChartification::ExpandCharts(CChartSet *charts)
{
   Heap<Halfedge *, double> heap;
   // clear chart map
   memset(chartMap, 0, sizeof(CChart *) * numFaces);
   // reset chart boundary
   boundary->Reset();
   // set last face to zero
   lastface = 0;

   // Initialize Heap
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // process single chart
      CChart *chart = *it;
      chart->Clear();

      // get seed information and create chart
      chartMap[chart->Seed()] = chart;

      // initialize boundary of chart
      Halfedge *loop = &halfedges[chart->Seed() * 3];

      // put all halfedges on heap
      for (int i=0; i<3; i++) {
         if (loop->sym) heap.push(loop, EdgeCost(loop));
         loop = loop->next;
      }
   }
   // local variables for chart growing
   Halfedge *h = NULL;
   double value = 0.0;

   // Charts-growing phase
   while(heap.pop(&h, &value)) {
      // continue with chart growing
      int f    = h->face;
      int fopp = h->sym->face;

      // check if we have a chart for opposite(f)
      if (chartMap[fopp] == NULL) {
         // add fopp to chart(F)
         chartMap[f]->Grow(faceNormals[fopp]);
         chartMap[fopp] = chartMap[f];
         
         // set last face value
         lastface = fopp;

         // remove E from chart_boundaries
         // for halfedges we have to erase both halfedge and symmetric halfedge
         boundary->Erase(h);

         // add halfedges of fopp belonging to boundary to chart-boundary
         Halfedge *loop = h->sym->next;
         for (int i=0; i<2; i++) {
            if (loop->sym && boundary->IsBoundary(loop)) heap.push(loop, EdgeCost(loop));
            loop = loop->next;
         }
         // repair heap
         //RepairHeap(chartMap[f], &heap);
      }
   }
}

// ----------------------------------------------------------------------------
// adds a face and boundary information to chart
// ----------------------------------------------------------------------------
void CChartification::AddInfoToChart(CChart *chart, Halfedge *h)
{
   // Create a face
   chart->AddFace(h->next->next->end, h->next->end, h->end);

   // check boundary info
   for (int i=0; i<3; i++) {
      // next halfedge
      h = h->next;
      // check if halfedge is boundary
      if (boundary->IsBoundary(h)) 
         chart->AddBoundary(h->end, h->next->next->end);
   }
}

// ----------------------------------------------------------------------------
// searches is vertex is already allocated
// ----------------------------------------------------------------------------
int CChartification::SearchVertexId(Halfedge *h, bool *visited)
{
   // Search first direction
   Halfedge *loop = h->next->sym;

   while (loop && loop!=h) {
      if (visited[loop->face]) return loop->end;
      loop = loop->next->sym;
   }

   // Search second direction
   if (h->sym) loop = h->sym->next->next;
   else        loop = NULL;

   while (loop && loop!=h) {
      if (visited[loop->face]) return loop->end;
      // next face
      if (loop->sym) loop = loop->sym->next->next;
      else           loop = NULL;
   }
   // nothing found
   return -1;
}

// ----------------------------------------------------------------------------
// creates the vertex info structure
// ----------------------------------------------------------------------------
void CChartification::CreateVertexStructure(CChartSet *charts, CChartBoundary *boundary)
{
   // heap for search algorithm
   Heap<Halfedge *, double> heap;
   Halfedge                 *loop;
   double                   value;

   // initialize markers
   bool *visited = new bool[numFaces];
   for (int i=0; i<numFaces; i++) visited[i] = false;

   // cut the halfedge mesh representation
   CutBoundaryPath(boundary);

   // initialize heap
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // process chart
      CChart *chart = *it;
      // get halfedge
      Halfedge *h = &halfedges[chart->Seed() * 3];

      // mark face as visited
      visited[h->face] = true;

      // process three halfedges
      for (int i=0; i<3; i++) {
         // next halfedge
         h = h->next;
         // map to a vertex and add it to heap
         h->end = chart->AddVertex(h->end);
         heap.push(h, 0.0);
      }
      // add information to chart
      AddInfoToChart(chart, h);
   }

   // apply a search over the dual graph of the mesh
   while (heap.pop(&loop, &value)) {
      // step to symmetric halfedge
      loop = loop->sym;

      // check if we've already been there
      if (!loop || visited[loop->face] || !chartMap[loop->face])
         continue;

      // get access to chart
      CChart *chart = chartMap[loop->face];

      // mark place as visited
      visited[loop->face] = true;

      // take to endpoints from symmetric halfedge
      loop->end = loop->sym->next->next->end;
      loop->next->next->end = loop->sym->end;

      // need a third endpoint
      loop = loop->next;

      // an beiden restlichen angrenzenden faces nach 
      // einem passenden Endpunkt suchen
      if (loop->sym && visited[loop->sym->face]) {
         // punkt an zweitem angrenzenden face gefunden
         loop->end = loop->sym->next->next->end;
         if (loop->next->next->end != loop->sym->end)
            printf ("heavy error occurred\n");
      }
      else if (loop->next->sym && visited[loop->next->sym->face]) {
         // punkt an dritten angrenzenden face gefunden
         loop->end = loop->next->sym->end;
      }
      else {
         // search for allocated vertex
         int index = SearchVertexId(loop, visited);
         // create a new one
         if (index < 0) loop->end = chart->AddVertex(loop->end);
         else           loop->end = index;
      }

      // add information to chart
      AddInfoToChart(chart, loop);

      // Kandidaten auf heap
      heap.push(loop, value + 1.0);
      heap.push(loop->next, value + 1.0);
   }
   delete[] visited;
}


// ----------------------------------------------------------------------------
// cuts along the boundary path, removing symmetrie info
// from boundary halfedges
// ----------------------------------------------------------------------------
void CChartification::CutBoundaryPath(CChartBoundary *boundary)
{
   for (int i=0; i<numHalfedges; i++) {
      // get halfedge
      Halfedge *h = &halfedges[i];
      // check for boundary
      if (h->sym && boundary->IsBoundary(h)) {
         // cut symmetrie info
         h->sym->sym = NULL;
         h->sym      = NULL;
      }
   }
}

// ----------------------------------------------------------------------------
// straighten the boundary
// ----------------------------------------------------------------------------
bool CChartification::Straighten(CChartBoundary *boundary)
{
   // über alle Flächen iterieren, falls eine Fläche 2 Boundary
   // edges besitzt, werden diese dem genüberliegenden Chart zugewiesen
   bool changed = false;

   for (int i=0; i<numFaces; i++) {
      // erste halbecke
      Halfedge *h     = &halfedges[i * 3];
      CChart   *chart = chartMap[i];

      if (!chart || i==chart->Seed())
         continue;

      for (int j=0; j<3; j++) {
         // Naechste Halbecke
         h = h->next;

         // Bedingung pruefen
         if (!(boundary->IsBoundary(h) && boundary->IsBoundary(h->next)))
            continue;
         
         // Symmetrische Kante vorhanden ?
         if (!(h->sym && h->next->sym))
            continue;
         
         CChart *ochart1 = chartMap[h->sym->face];
         CChart *ochart2 = chartMap[h->next->sym->face];

         if (ochart1 && ochart2 && ochart1 == ochart2 && chart!=ochart1) {
            // Error Faktor bestimmen
            Vec3 &nrm1 = faceNormals[h->face];
            Vec3 &nrm2 = faceNormals[h->sym->face];
            
            // get Factor for difference in normals
            double err  = 1.0 - (nrm1 * nrm2) / sqrt(nrm1.length2() * nrm2.length2());

            if (err < 0.3)
            {
               // korrigieren
               chartMap[h->face] = ochart1;
               // boundary struktur anpassen
               boundary->SetBoundaryFlag(h->next->next, true);
               boundary->Erase(h);

               changed = true;

               // Abbruch
               break;
            }
         }
      }
   }
   return changed;
}

// ----------------------------------------------------------------------------
// allocates enough memory for the operations
// ----------------------------------------------------------------------------
void CChartification::Allocate(CChartSet  *charts,
                                    int   *numPoints,
                                    Vec3 **points,
                                    rgbf **rgb_values,
                                    int   *numFaces,
                                    Face **faces,
                                    int   *numEdges,
                                    Edge **edges
                              )
{
   // allocate space
   *numPoints = 0;
   *numFaces  = 0;
   *numEdges  = 0;

   // get overall count of elements
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // add chart sizes
      *numPoints += (*it)->VertexCount();
      *numFaces  += (*it)->FaceCount();
      *numEdges  += (*it)->BoundaryCount();
   }

   // allocate it
   *points     = new Vec3[*numPoints];
   *rgb_values = new rgbf[*numPoints];
   *faces      = new Face[*numFaces];
   *edges      = new Edge[*numEdges];
}

// ----------------------------------------------------------------------------
// sets a color value for each different chart
// ----------------------------------------------------------------------------
void CChartification::Colorize(CChartSet *charts, bool tex,
                                   int   *numPoints,
                                   Vec3 **points,
                                   rgbf **rgb_values,
                                   int   *numFaces,
                                   Face **faces,
                                   int   *numEdges,
                                   Edge **edges
                              )
{
   // allocate space
   Allocate(charts, numPoints, points, rgb_values, numFaces, faces, numEdges, edges);

   // initialize colorfactory
   CColorFactory::Initialize();

   // offsets for combining the charts
   int pntoff = 0;
   int facoff = 0;
   int edgoff = 0;

   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // process chart
      CChart *chart = *it;

      // generate a html compatible color
      rgbf col = CColorFactory::HtmlColor();

      // copy and map from vertex to vec3
      Vertex *cvertices = chart->Vertices();
      //
      for (int i=0; i<chart->VertexCount(); i++) {
         // fill reference
         Vec3 &pnt   = (*points)[i + pntoff];
         Vertex &vtx = cvertices[i];
         // map from vertex to vec3 index
         if (!tex) pnt = this->points[vtx.index];
         else {
            // write texture als coord
            pnt.x = vtx.u;
            pnt.y = vtx.v;
            pnt.z = 0.0f;
         }
         // write color for point
         (*rgb_values)[i + pntoff] = col;
      }

      // copy and remap faces
      Face *cfaces = chart->Faces();
      //
      for (int i=0; i<chart->FaceCount(); i++) {
         // copy face
         (*faces)[i + facoff].a = cfaces[i].a + pntoff;
         (*faces)[i + facoff].b = cfaces[i].b + pntoff;
         (*faces)[i + facoff].c = cfaces[i].c + pntoff;         
      }

      // copy and remap edges
      Edge *cedges = chart->Boundaries();
      //
      for (int i=0; i<chart->BoundaryCount(); i++) {
         // copy edge
         (*edges)[i + edgoff].a = cedges[i].a + pntoff;
         (*edges)[i + edgoff].b = cedges[i].b + pntoff;
      }

      // increase offset
      pntoff = pntoff + chart->VertexCount();
      facoff = facoff + chart->FaceCount();
      edgoff = edgoff + chart->BoundaryCount();
   }
}


// ----------------------------------------------------------------------------
// creates a chart visualization
// ----------------------------------------------------------------------------
void CChartification::TextureVisualization(int numDesiredCharts,
                                    int   *numPoints,
                                    Vec3 **points,
                                    rgbf **rgb_values,
                                    int   *numFaces,
                                    Face **faces,
                                    int   *numEdges,
                                    Edge **edges
                                    )
{
   // get the charts
   CChartSet *charts = new CChartSet(this->points, this->faceNormals, this->halfedges);
   Chartify(numDesiredCharts, charts);

   CConformalMapper mapper(this->points);
   mapper.Parametrize(charts);
 
   printf("Atlas packing\n");
   CAtlasPacker packer(this->points, 256);
   packer.Pack(charts);

   Colorize(charts, true, numPoints, points, rgb_values, numFaces, faces, numEdges, edges);

   // delete chart set
   delete charts;
}

// ----------------------------------------------------------------------------
// creates a chart visualization
// ----------------------------------------------------------------------------
void CChartification::Visualization(int numDesiredCharts,
                                    int   *numPoints,
                                    Vec3 **points,
                                    rgbf **rgb_values,
                                    int   *numFaces,
                                    Face **faces,
                                    int   *numEdges,
                                    Edge **edges
                                    )
{
   // get the charts
   CChartSet *charts = new CChartSet(this->points, this->faceNormals, this->halfedges);
   Chartify(numDesiredCharts, charts);

   Colorize(charts, false, numPoints, points, rgb_values, numFaces, faces, numEdges, edges);

   // delete chart set
   delete charts;
}


// ----------------------------------------------------------------------------
// Gets Vertices from Chartset representation
// ----------------------------------------------------------------------------
void CChartification::GetVertices(CChartSet *charts, int *numVertices, Vertex **vertices)
{
   // initialize
   *numVertices  = 0;

   // get overall count of elements
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // add chart sizes
      *numVertices  += (*it)->VertexCount();
   }
   // allocate it
   *vertices = new Vertex[*numVertices];

   // offsets for combining the charts
   int pntoff = 0;

   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // process chart
      CChart *chart = *it;

      // copy and map from vertex to vec3
      Vertex *cvertices = chart->Vertices();
      //
      for (int i=0; i<chart->VertexCount(); i++) {
         // copy vertex
         (*vertices)[pntoff + i] = cvertices[i];
      }

      // increase offset
      pntoff = pntoff + chart->VertexCount();
   }
}


// ----------------------------------------------------------------------------
// Gets Faces from Chartset representation
// ----------------------------------------------------------------------------
void CChartification::GetFaces(CChartSet *charts, int *numFaces, Face **faces)
{
   // initialize
   *numFaces  = 0;

   // get overall count of elements
   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // add chart sizes
      *numFaces  += (*it)->FaceCount();
   }
   // allocate it
   *faces      = new Face[*numFaces];

   // offsets for combining the charts
   int pntoff = 0;
   int facoff = 0;


   for (CChartSet::iterator it=charts->Begin(); it!=charts->End(); it++) {
      // process chart
      CChart *chart = *it;

      // copy and remap faces
      Face *cfaces = chart->Faces();
      //
      for (int i=0; i<chart->FaceCount(); i++) {
         // copy face
         (*faces)[i + facoff].a = cfaces[i].a + pntoff;
         (*faces)[i + facoff].b = cfaces[i].b + pntoff;
         (*faces)[i + facoff].c = cfaces[i].c + pntoff;         
      }

      // increase offsets
      pntoff = pntoff + chart->VertexCount();
      facoff = facoff + chart->FaceCount();
   }
}


// ----------------------------------------------------------------------------
// creates a chart visualization
// ----------------------------------------------------------------------------
void CChartification::CreateVertices(// incoming parameters
                                     int numDesiredCharts,

                                     // outgoing parameters
                                     int     *outNumVertices,
                                     Vertex **outVertices,
                                     int     *outNumFaces,
                                     Face   **outFaces
                                    )
{
   // get the charts
   CChartSet *charts = new CChartSet(this->points, this->faceNormals, this->halfedges);
   Chartify(numDesiredCharts, charts);

   CConformalMapper mapper(this->points);
   mapper.Parametrize(charts);
 
   printf("Atlas packing\n");
   CAtlasPacker packer(this->points, 256);
   packer.Pack(charts);

   GetFaces(charts, outNumFaces, outFaces);
   GetVertices(charts, outNumVertices, outVertices);

   // delete chart set
   delete charts;
}
} // end namespace

