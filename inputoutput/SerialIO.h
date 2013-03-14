/***************************************************************************
  SerialIO.h  -  a class to write a file per processor 
  -------------------

developers           : Stefano Markidis, Giovanni Lapenta

 ***************************************************************************/

#ifndef SerialIO_H
#define SerialIO_H



#include <iostream>
#include <sstream>
#include <string>
#include <fstream>


using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;


/**
 * Library to write a file per processor
 *
 */
/** write to an ASCII file Electrostatic Potential */
inline void writePHIascii1D(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  int nxc = grid->getNXC();

  for (int i = 1; i < grid->getNXC() - 1; i++) {
    my_file << grid->getXC(i, 0, 0) << "\t " << (field->getPHI())[i][0][0] << endl;
  }
  my_file.close();


}
inline void writeRHOascii1D(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  int nxc = grid->getNXC();

  for (int i = 0; i < grid->getNXN(); i++) {
    my_file << grid->getXN(i, 0, 0) << "\t " << field->getRHOn(i, 0, 0) << endl;
  }
  my_file.close();
}
inline void writeRHOascii2D(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  int nxc = grid->getNXC();

  for (int i = 1; i < grid->getNXN() - 1; i++) {
    for (int j = 1; j < grid->getNYN() - 1; j++) {
      my_file << field->getRHOn(i, j, 0) << "\t";
    }
    my_file << endl;
  }
  my_file.close();
}
inline void writeRHOascii2Del(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  int nxc = grid->getNXC();

  for (int i = 1; i < grid->getNXN() - 1; i++) {
    for (int j = 1; j < grid->getNYN() - 1; j++) {
      my_file << field->getRHOns(i, j, 0, 0) << "\t";
    }
    my_file << endl;
  }
  my_file.close();
}
inline void writeRHOascii2Dion(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  int nxc = grid->getNXC();

  for (int i = 1; i < grid->getNXN() - 1; i++) {
    for (int j = 1; j < grid->getNYN() - 1; j++) {
      my_file << field->getRHOns(i, j, 0, 1) << "\t";
    }
    my_file << endl;
  }
  my_file.close();
}
/** write an ASCII with Ex */
inline void writeExascii1D(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./resultsOsc/" + filename + "_cycle" + cc.str() + "_proc" + ss.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  for (int i = 0; i < (grid->getNXN()); i++) {
    my_file << field->getRHOns(i, 0, 0, 0) << endl;
  }
  my_file.close();
}
/** write an ASCII with the heat flux */
inline void writeJx1D(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  for (int i = 0; i < (grid->getNXN()); i++) {
    my_file << field->getJx(i, 0, 0) << endl;
  }
  my_file.close();
}
/** write an ASCII with the heat flux */
inline void writeHeatFluxascii1D(string filename, int myrank, int cycle, Grid * grid, Field * field) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./resultsOsc/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  for (int i = 0; i < (grid->getNXN()); i++) {
    my_file << field->getJy(i, 0, 0) << endl;
  }
  my_file.close();
}
/** write to an ASCII file  particles position */
inline void writeParticles1D(string filename, int myrank, int cycle, Particles * part) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./results/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  for (int i = 0; i < part->getNOP(); i++) {
    my_file << part->getX(i) << "\t" << part->getU(i) << endl;
  }
  my_file.close();
}
/** write to an ASCII file  particles position */
inline void writeParticles1D(string filename, int myrank, int cycle, Particles * part, int p) {
  // stream file to be opened and managed
  string temp;
  stringstream ss;
  stringstream cc;
  ss << myrank;
  cc << cycle;
  temp = "./resultsOsc/" + filename + "_proc" + ss.str() + "_cycle" + cc.str();
  temp += ".txt";
  cout << "Opening file: " << temp << endl;
  ofstream my_file(temp.c_str());
  for (int i = 0; i < part->getNOP(); i = i + p) {
    my_file << part->getX(i) << "\t" << part->getU(i) << endl;
  }
  my_file.close();
}


#endif
