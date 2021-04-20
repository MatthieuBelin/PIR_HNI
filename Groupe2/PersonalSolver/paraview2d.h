// check dimension
#if dimension != 2
#error ERROR: dimension != 2 
#endif

// includes for mkdir and time
#include <sys/stat.h>
#include <time.h>

// result's folder 
static char dirname[80], pvd[80];

// init result's folder
event defaults (i = 0){

// get dirname
  time_t rawtime;  time( &rawtime );
  struct tm *info; info = localtime( &rawtime );
  strftime(dirname,80,"%Y%m%d-%H%M%S", info);
// init pvd file
  sprintf(pvd, "%s/result.pvd", dirname);

#if _MPI
  if (pid()) // !master
    return 0;
#endif
  
//checked if directory is exists
  struct stat st;
  if (stat(dirname, &st) == -1)
    mkdir(dirname, 0700);
  
  FILE * fp = fopen(pvd, "w"); assert(fp);
  fputs ("<?xml version=\"1.0\"?>\n"
	 "<VTKFile type=\"Collection\" version=\"1.0\" \
byte_order=\"LittleEndian\">\n", fp);
  fputs ("\t <Collection>\n", fp);
  fclose(fp);
}

// finalize result's folder
event finalise(i = end){

#if _MPI
  if (pid()) // !master
    return 0;
#endif
  
  // open the file result.pvd and write 
  FILE * fp = fopen(pvd, "a"); assert(fp);
  fputs ("\t </Collection>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fclose(fp);
}

// prototype
void output_pvtu (scalar * list, vector * vlist, FILE * fp, char * vtufiles);
void output_vtu (scalar * list, vector * vlist, FILE * fp);

struct OutParaview {
  scalar * slist;
  vector * vlist;
  bool listing;
};

void output_paraview(struct OutParaview p) {
  // default value for optional arguments
  scalar * slist = p.slist;
  vector * vlist = p.vlist;
  bool listing = p.listing;
  //
  double time = t;
  static int nf = -1; // number of files
  char vtu[80], name[80]; // file name to open
  FILE * fp;     // flush
  
  nf++; // new file
  // define vtu file name
#if _MPI
  sprintf(vtu, "ts_n%4.4d_p%4.4d.vtu", nf, pid());
#else
  sprintf(vtu, "ts_n%4.4d.vtu", nf);
#endif
  // write vtu file
  sprintf(name, "%s/%s", dirname, vtu);
  fp = fopen(name, "w"); assert(fp);
  output_vtu (slist, vlist, fp);
  fclose (fp);

#if _MPI
  // master proc writes pvtu and pvd file
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid()) // !master
    return; 
  
  char pvtu[80];
  // redefine vtu file for pvtu
  sprintf(vtu, "ts_n%4.4d", nf); 
  sprintf(pvtu, "%s.pvtu", vtu);
  sprintf(name, "%s/%s", dirname, pvtu);
  fp = fopen(name, "w"); assert(fp);
  output_pvtu (slist, vlist, fp, vtu);
  fclose (fp);
  // redefine vtk file for pvd
  sprintf(vtu, "%s", pvtu);
#endif 
  
  // add vtu file into pvd file
  if (listing)
    fprintf(stdout, "# t = %f \n", time);
  fp = fopen(pvd, "a"); assert(fp);
  fprintf (fp,"\t\t <DataSet timestep=\"%f\" file=\"%s\"/>\n", time, vtu);
  fclose(fp);
}

void output_pvtu (scalar * slist, vector * vlist, FILE * fp, char * vtufiles) {
  fputs ("<?xml version=\"1.0\"?>\n"
	 "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" \
byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
  fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
  
  if (slist)
    for (scalar s in slist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" \
format=\"appended\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }

  if (vlist)
    for (vector v in vlist) {
      char vname[32]; strcpy(vname, v.x.name); vname[strlen(vname) - 2] = '\0'; 
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" \
NumberOfComponents=\"3\" Name=\"%s\" format=\"appended\">\n", vname);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }

  fputs ("\t\t\t </PCellData>\n", fp);
  fputs ("\t\t\t <PPoints>\n", fp);
  fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" \
format=\"ascii\">\n", fp);
  fputs ("\t\t\t\t </PDataArray>\n", fp);
  fputs ("\t\t\t </PPoints>\n", fp);

  for (int i = 0; i < npe(); i++)
    fprintf (fp, "<Piece Source=\"%s_p%4.4d.vtu\"/> \n", vtufiles, i);

  fputs ("\t </PUnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
}

/*
  This function writes one XML VTK file per PID process of type unstructured grid
  (*.vtu) which can be read using Paraview. File stores scalar and vector fields
  defined at the center points. Results are recorded on binary format. If one writes
  one *.vtu file per PID process this function may be combined with
  output_pvtu() above to read in parallel. Tested in (quad- and oct-)trees
  using MPI. Also works with solids (when not using MPI).
*/
void output_vtu (scalar * slist, vector * vlist, FILE * fp) {

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
      //marker[] = _k;
    marker[] = no_points;
    no_points += 1;
  }
  foreach(){
    no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
	 "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" \
byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
	   no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  int count = 0;
    
  if (slist)
    for (scalar s in slist) {
      fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" \
format=\"appended\" offset=\"%d\">\n", s.name,count);
      count += ((no_cells)+1)*8;
      fputs ("\t\t\t\t </DataArray>\n", fp);
    }

  if (vlist)
    for (vector v in vlist) {
      char vname[32]; strcpy(vname, v.x.name); vname[strlen(vname) - 2] = '\0'; 
      fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" \
NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", vname,count);
      count += ((no_cells*3)+1)*8;
      fputs ("\t\t\t\t </DataArray>\n", fp);
    }

  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" \
format=\"appended\" offset=\"%d\">\n",count);
  count += ((no_points*3)+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" \
format=\"ascii\">\n", fp);
  foreach(){
    fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], 
	     marker[0,1]);
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" \
format=\"ascii\">\n", fp);
  for (int i = 1; i < no_cells+1; i++){
    fprintf (fp, "%d \n", i*4);
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" \
format=\"ascii\">\n", fp);
  foreach(){
    fputs ("9 \n", fp);
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);

  unsigned long long block_len;

  if (slist) {
    block_len=no_cells*8;
    for (scalar s in slist) {
      fwrite (&block_len, sizeof (unsigned long long), 1, fp);
      foreach()
	fwrite (&val(s), sizeof (double), 1, fp);
    }
  }
      
  if (vlist) {
    block_len=no_cells*8*3;
    for (vector v in vlist) {
      fwrite (&block_len, sizeof (unsigned long long), 1, fp);
      foreach(){
	fwrite (&val(v.x), sizeof (double), 1, fp);
	fwrite (&val(v.y), sizeof (double), 1, fp);
	double vz=0;
	fwrite (&vz, sizeof (double), 1, fp);
      }
    }
  }
  
  block_len=no_points*8*3;
  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
  foreach_vertex(){
    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
}

