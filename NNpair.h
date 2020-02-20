#ifndef NNPAIR_H
#define NNPAIR_H

#include <string>

class NNpair
{
private:
  std::string q_id; //to id tou q
  std::string p_id; //to id tou kontinoterou geitona sto dataset
  double distance;

public:
  NNpair(){};
  NNpair(std::string q, std::string p); // o conustructor gia parametropoihsh
  std::string getq_id();
  int getq_id_as_int();
  void setq_id(std::string idd);
  std::string getp_id();
  int getp_id_as_int();
  void setp_id(std::string idd);
  double get_distance();
  void set_distance(double dis);
};

#endif
