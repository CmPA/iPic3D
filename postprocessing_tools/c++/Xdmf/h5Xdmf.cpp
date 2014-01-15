#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

class c_time{
  public:
    c_time();
    ~c_time() {;} ;
    void initcycledata(int b, int e, int s, int n);
    int  getncycle(){return ncycle;};
    int  geticycle(int i){return (i*stp);};
    int  getbeg(){return beg;};
    int  getend(){return end;};
    int  getstp(){return stp;};

  private:
    int beg;
    int end;
    int stp;
    int ncycle;
};

c_time::c_time(){
  beg = 0;
  end = 0;
  stp = 0;
  ncycle = 0;
}

void c_time::initcycledata(int b, int e, int s, int n){
  beg    = b;
  end    = e;
  stp    = s;
  ncycle = n;
}

class c_grid{
  public:
    c_grid();
    ~c_grid() {;};
    void setdimorigin(int nd, double *org);
    void setdelta(double *delta);

    double getOx(){return Ox;};
    double getOy(){return Oy;};
    double getOz(){return Oz;};

    double getdx(){return dx;};
    double getdy(){return dy;};
    double getdz(){return dz;};

  private:
    int    ndim;
    double dx;
    double dy;
    double dz;
    double Ox;
    double Oy;
    double Oz;
};

c_grid::c_grid(){
  ndim = 2;
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;
  Ox = 0.0;
  Oy = 0.0;
  Oz = 0.0;
}

void c_grid::setdimorigin(int nd, double *org){
  ndim = nd;

  Ox = org[0];
  Oy = org[1];
  if (ndim == 3) Oz = org[2];

}

void c_grid::setdelta(double *delta){

  dx = delta[0];
  dy = delta[1];
  if (ndim == 3) dz = delta[2];

}

class c_flds{

  public:
    c_flds() {;};
    ~c_flds();
    void initfields(int nf);
    void setfieldname(int i, string dtset, string fieldname);
    int  getnfields(){return nfields;};
    string getidtset (int i){return dataset[i];};
    string getidtname(int i){return dtname [i];};
  private:
    int     nfields;
    string *dataset;
    string *dtname;
};

void c_flds::initfields(int nf){
  nfields = nf;
  dataset = new string[nf];
  dtname  = new string[nf];
}

void c_flds::setfieldname(int i, string dtset, string fieldname){
  dataset[i] = dtset;
  dtname [i] = fieldname;
}

c_flds::~c_flds(){
  delete [] dataset;
  delete [] dtname;
}

class xdmfile{

  public:
    xdmfile();
    ~xdmfile() {;};
    void writexdmf();
    void readinput(string infile);
    void setgeometry(void);
    void setcycles(void);
    void setfields(void);
    string h5filename;

    int get_cbeg(){return cbeg;};
    int get_cend(){return cend;};
    int get_cstp(){return cstp;};

    int set_cbeg(int i){cbeg = i;};
    int set_cend(int i){cend = i;};
    int set_cstp(int i){cstp = i;};

  private:

    c_time time;
    c_grid grid;
    c_flds fields;

    int nfields;
    int ndim;
    int ncx;
    int ncy;
    int ncz;

    int cbeg;
    int cend;
    int cstp;

    double Lx;
    double Ly;
    double Lz;
  
};

xdmfile::xdmfile(){
  int nfields = 1;
  int ndim    = 2;
  int ncx     = 0;
  int ncy     = 0;
  int ncz     = 0;
  int cbeg    = 0;
  int cend    = 0;
  int cstp    = 0;

  double Lx   = 0.0;
  double Ly   = 0.0;
  double Lz   = 0.0;
}

void xdmfile::setgeometry(){

  double *arr;
  arr = new double[ndim];

  arr[0] = 0.0;
  arr[1] = 0.0;
  if (ndim == 3) arr[2] = 0.0;

  grid.setdimorigin(ndim,arr);

  arr[0] = Lx/double(ncx);
  arr[1] = Ly/double(ncy);
  if (ndim == 3) arr[2] = Lz/double(ncz);

  grid.setdelta(arr);

  delete [] arr;

}

void xdmfile::setcycles(){
  time.initcycledata(get_cbeg(), get_cend(), get_cstp(), int((get_cend()-get_cbeg())/get_cstp()) + 1);
}

void xdmfile::setfields(){
}

void xdmfile::readinput(string infile){

  ifstream     inp(infile.c_str());

  int    beg, end, stp, nf;
  string fieldline;
  string dtset;
  string fieldname;

  inp >> h5filename;
  inp >> ndim;
  inp >> beg;
  inp >> end;
  inp >> stp;
  inp >> Lx;
  inp >> Ly;
  inp >> Lz;
  inp >> ncx;
  inp >> ncy;
  inp >> ncz;
  inp >> nf;

  fields.initfields(nf);

  for (int i = 0; i< fields.getnfields(); i++){
    inp >> dtset;
    inp >> fieldname;
    fields.setfieldname(i, dtset, fieldname);
  }

  set_cbeg(beg);
  set_cend(end);
  set_cstp(stp);

  if (ndim < 2 || ndim > 3) abort();

  setgeometry();
  setcycles();
  setfields();

}

void xdmfile::writexdmf(){

  stringstream ss;

  ss << h5filename << ".xmf";

  ofstream xdmfout;
  xdmfout.open(ss.str().c_str());

  // ------
  // Header
  // ------

  xdmfout << "<?xml version=\"1.0\" ?>" << endl;
  xdmfout << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"> " << endl;
  xdmfout << "<Xdmf Version=\"2.0\"> " << endl;
  xdmfout << "  <Domain> " << endl;

  // ------------------------------
  // Write time sequence parameters
  // ------------------------------

  xdmfout << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\"> " << endl;
  xdmfout << " " << endl;
  xdmfout << "      <Time TimeType=\"List\"> " << endl;

  xdmfout << "        <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"" << time.getncycle() << "\"> " << endl;
  xdmfout << "          ";

  for (int i = 0; i < time.getncycle(); i++)
    xdmfout << time.geticycle(i) << " ";

  xdmfout << endl;

  xdmfout << "        </DataItem> " << endl;
  xdmfout << "      </Time> " << endl;

  // ----------------
  // Loop over all the cycles
  // ----------------

  for (int t = time.getbeg(); t <= time.getend(); t+=time.getstp()){

    xdmfout << endl;

    stringstream filenmbr;
    string       filename;

    filenmbr << setfill('0') << setw(5) << t;
    filename = h5filename + "_" + filenmbr.str() + ".h5";

    if (ndim == 2){

      xdmfout << "      <Grid Name=\"Structured mesh\" GridType=\"Uniform\"> " << endl;
      xdmfout << "        <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << ncy+1 << " " << ncx+1 << "\"/> " << endl;
      xdmfout << "        <Geometry GeometryType=\"ORIGIN_DXDY\"> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"2\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << setprecision(2) << fixed << grid.getOy() << " " << setprecision(2) << fixed << grid.getOx() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"2\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << grid.getdy() << " " << grid.getdx() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "        </Geometry> " << endl;

    // ----------------
    // Write each field
    // ----------------

      for (int i = 0; i < fields.getnfields(); i++){

        xdmfout << "        <Attribute Name=\"" << fields.getidtname(i) << "\" Center=\"Cell\"> " << endl;
        xdmfout << "          <DataItem Format=\"HDF\" Dimensions=\"" << ncy << " " << ncx << "\" NumberType=\"Float\"> " << endl;
        xdmfout << "          " <<  filename << ":" << fields.getidtset(i) << endl;
        xdmfout << "          </DataItem> " << endl;
        xdmfout << "        </Attribute> " << endl;

      }
    }
    else{

      xdmfout << "      <Grid Name=\"Structured mesh\" GridType=\"Uniform\"> " << endl;
      xdmfout << "        <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << ncx+1 << " " << ncy+1 << " " << ncz+1 << "\"/> " << endl;
      xdmfout << "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\"> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"3\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << setprecision(2) << fixed << grid.getOx() << " " << setprecision(2) << fixed << grid.getOy() << " " << setprecision(2) << fixed << grid.getOz() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"3\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << grid.getdx() << " " << grid.getdy() << " " << grid.getdz() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "        </Geometry> " << endl;

    // ----------------
    // Write each field
    // ----------------

      for (int i = 0; i < fields.getnfields(); i++){

        xdmfout << "        <Attribute Name=\"" << fields.getidtname(i) << "\" Center=\"Cell\"> " << endl;
        xdmfout << "          <DataItem Format=\"HDF\" Dimensions=\"" << ncx << " " << ncy << " " << ncz << "\" NumberType=\"Float\"> " << endl;
        xdmfout << "          " <<  filename << ":" << fields.getidtset(i) << endl;
        xdmfout << "          </DataItem> " << endl;
        xdmfout << "        </Attribute> " << endl;

      }

    }

    xdmfout << "      </Grid> " << endl;

  }  // End cycles loop

  // ----------------------
  // Close Time grid/Domain/Xdmf tags
  // ----------------------

  xdmfout << endl;
  xdmfout << "    </Grid> " << endl;
  xdmfout << "  </Domain> " << endl;
  xdmfout << "</Xdmf> " << endl;

  // ----------
  // Close file
  // ----------

  xdmfout.close();

}

int main(void) {

  xdmfile xdmf;
  string infile = "h5Xdmf.inp";
  xdmf.readinput(infile);
  xdmf.writexdmf();

  return 0;

}
