#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>

#include "ConfigFile.h"

using namespace std;

struct filefield{
  string fieldname;
  string dtset;

  filefield & operator=(const filefield & orig) {
    fieldname = orig.fieldname;
    dtset     = orig.dtset;
  }
};

inline ostream & operator<<(ostream & os, const filefield & t){
  os << t.fieldname << " " << t.dtset;
  return os;
}

inline istream & operator>>(istream & is, filefield & t){
  is >> t.fieldname >> t.dtset;
  return is;
}

class c_time{
  public:
    c_time();
    ~c_time() {;} ;
    void initcycledata(int b, int e, int s, int n);
    int  getncycle(){return ncycle;};
    int  geticycle(int i){return (beg+i*stp);};
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
    void setlength(double *L);

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
    double Lx;
    double Ly;
    double Lz;
};

c_grid::c_grid(){
  ndim = 2;
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;
  Ox = 0.0;
  Oy = 0.0;
  Oz = 0.0;
  Lx = 0.0;
  Ly = 0.0;
  Lz = 0.0;
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

void c_grid::setlength(double *L){

  Lx = L[0];
  Ly = L[1];
  if (ndim == 3) Lz = L[2];

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
    void readinput_org(string infile);
    void setgeometry(void){};
    void setgeometry(int, double*, double*);
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
}

void xdmfile::setgeometry(int ndim, double *O, double *L){

  double *arr;
  arr = new double[ndim];

  arr[0] = L[0]/double(ncx);
  arr[1] = L[1]/double(ncy);
  if (ndim == 3) arr[2] = L[2]/double(ncz);

  grid.setdimorigin(ndim, O);
  grid.setlength(L);
  grid.setdelta(arr);

  delete [] arr;

}

void xdmfile::setcycles(){
  time.initcycledata(get_cbeg(), get_cend(), get_cstp(), int((get_cend()-get_cbeg())/get_cstp()) + 1);
}

void xdmfile::setfields(){
}

void xdmfile::readinput(string infile){
  ConfigFile inp(infile);
  int        beg, end, stp, nf;
  double     Ox, Oy, Oz;
  double     Lx, Ly, Lz;

  h5filename = inp.read<string>("H5FileName");
  ndim       = inp.read<int>   ("NumberOfDimensions");
  beg        = inp.read<int>   ("FirstCycle");
  end        = inp.read<int>   ("LastCycle");
  stp        = inp.read<int>   ("StepBetweenCycles");
  Ox         = inp.read<double>("OriginInX");
  Oy         = inp.read<double>("OriginInY");
  Oz         = inp.read<double>("OriginInZ");
  Lx         = inp.read<double>("LenghtInX");
  Ly         = inp.read<double>("LenghtInY");
  Lz         = inp.read<double>("LenghtInZ");
  ncx        = inp.read<int>   ("NumberOfCellsInX");
  ncy        = inp.read<int>   ("NumberOfCellsInY");
  ncz        = inp.read<int>   ("NumberOfCellsInZ");
  nf         = inp.read<int>   ("NumberOfFields");

  fields.initfields(nf);

  for (int i = 0; i< fields.getnfields(); i++){
    filefield ff;
    stringstream ss;
    ss << i+1;
    string s_i = ss.str();
    ff = inp.read<filefield>("Field_"+s_i);
    fields.setfieldname(i, ff.dtset, ff.fieldname);
  }

  set_cbeg(beg);
  set_cend(end);
  set_cstp(stp);

  if (ndim < 2 || ndim > 3) abort();

  double *O = new double [ndim];
  double *L = new double [ndim];

  O[0] = Ox;
  O[1] = Oy;
  if (ndim==3) O[2] = Oz;
  L[0] = Lx;
  L[1] = Ly;
  if (ndim==3) L[2] = Lz;

  setgeometry(ndim, O, L);
  setcycles();
  setfields();

  delete [] O;
  delete [] L;

}

void xdmfile::readinput_org(string infile){

  ifstream     inp(infile.c_str());
  double     Ox, Oy, Oz;
  double     Lx, Ly, Lz;

  int    beg, end, stp, nf;
  string fieldline;
  string dtset;
  string fieldname;

  inp >> h5filename;
  inp >> ndim;
  inp >> beg;
  inp >> end;
  inp >> stp;
  inp >> Ox;
  inp >> Oy;
  inp >> Oz;
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

    filenmbr << setfill('0') << setw(6) << t;
    filename = h5filename + "_" + filenmbr.str() + ".h5";

    if (ndim == 2){

      xdmfout << "      <Grid Name=\"Structured mesh\" GridType=\"Uniform\"> " << endl;
      xdmfout << "        <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << ncy+1 << " " << ncx+1 << "\"/> " << endl;
      xdmfout << "        <Geometry GeometryType=\"ORIGIN_DXDY\"> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"2\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << setprecision(4) << fixed << grid.getOy() << " " << setprecision(4) << fixed << grid.getOx() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"2\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << grid.getdy() << " " << grid.getdx() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "        </Geometry> " << endl;

    // ----------------
    // Write each field
    // ----------------

      for (int i = 0; i < fields.getnfields(); i++){

        xdmfout << "        <Attribute Name=\"" << fields.getidtname(i) << "\" Center=\"Node\"> " << endl;
        xdmfout << "          <DataItem Format=\"HDF\" Dimensions=\"" << ncy+1 << " " << ncx+1 << "\" NumberType=\"Float\"> " << endl;
        xdmfout << "          " <<  filename << ":" << fields.getidtset(i) << endl;
        xdmfout << "          </DataItem> " << endl;
        xdmfout << "        </Attribute> " << endl;

      }
    }
    else{

      xdmfout << "      <Grid Name=\"Structured mesh\" GridType=\"Uniform\"> " << endl;
      xdmfout << "        <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << ncz+1 << " " << ncy+1 << " " << ncx+1 << "\"/> " << endl;
      xdmfout << "        <Geometry GeometryType=\"ORIGIN_DXDYDZ\"> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"3\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << setprecision(4) << fixed << grid.getOz() << " " << setprecision(4) << fixed << grid.getOy() << " " << setprecision(4) << fixed << grid.getOx() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "          <DataItem Format=\"XML\" Dimensions=\"3\" NumberType=\"Float\"> " << endl;
      xdmfout << "          " << grid.getdz() << " " << grid.getdy() << " " << grid.getdx() << endl;
      xdmfout << "          </DataItem> " << endl;
      xdmfout << "        </Geometry> " << endl;

    // ----------------
    // Write each field
    // ----------------

      for (int i = 0; i < fields.getnfields(); i++){

        xdmfout << "        <Attribute Name=\"" << fields.getidtname(i) << "\" Center=\"Node\"> " << endl;
        xdmfout << "          <DataItem Format=\"HDF\" Dimensions=\"" << ncz+1 << " " << ncy+1 << " " << ncx+1 << "\" NumberType=\"Float\"> " << endl;
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

int main(int argc, char *argv[]) {

  if (argc==2) {
    xdmfile xdmf;
    string infile(argv[1]);
    xdmf.readinput(infile);
    xdmf.writexdmf();
  }
  else {
    cout << " USE: h5Xdmf <inputfile>" << endl;
  }

  return 0;

}
