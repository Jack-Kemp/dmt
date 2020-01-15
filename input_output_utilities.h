#ifndef INPUT_OUTPUT_UTILITIES_H
#define INPUT_OUTPUT_UTILITIES_H

#include<iostream>
#include<fstream>
#include<string>
#include<map>
#include<functional>
#include<vector>

struct Flags {
  bool periodic, antiperiodic, odd_sector, eentropy_growth;
  bool excited_state, correlations_decay, eentropy_c, plaqexp;
  bool leg_inversion; bool rung_basis;
};

template <class A, class B>
std::string construct_label(const std::string& base,
                            const std::string& vary,
                            Flags f, A params, B params_int){
 std::string label = base;
  if (f.periodic) label +=  "_periodic";
  else if (f.antiperiodic) label +=  "_antiperiodic";
  else label +=  "_open";
  if (not f.rung_basis){
    if (not f.leg_inversion){
      if (f.odd_sector) label += "_odd";
      else label += "_even";
    }
  }
  for (auto const& x : params)
    {
      label += "_" + x.first + "_";
      label += (x.first == vary ? "vary" : std::to_string(x.second));
    }
  for (auto const& x : params_int)
    {
      label += "_" + x.first + "_";
      label += (x.first == vary ? "vary" : std::to_string(x.second));
    }
  return label;
}

template<typename X>
void readfromFile (std::string iName, X& x){
        std::ifstream file;
        file.open(iName.c_str());
        if (file)
        	for  (int i =0; i<x.size()&& file >> x[i]; i++);
        else std::cout<<"Could not open file!"<<std::endl;
}

template<typename X>
void writeToFile (std::string iName,const  X& x, bool rewrite){
        std::ofstream outfile;
        if (rewrite)
          outfile.open(iName.c_str(), std::ios::trunc);
        else
          outfile.open(iName.c_str(), std::ios::app);
        if (outfile){
        outfile << std::string("# ") + iName<<'\n';
        for (int i =0; i<x.size();i++)
            outfile<< x[i]<<'\n';
        outfile.flush();
        outfile.close();
        }
}

template<typename X, typename Y>
void writeToFile (std::string iName,const  X& x, const  Y& y, bool rewrite){
        std::ofstream outfile;
        if (x.size()==y.size()){
        if (rewrite)
          outfile.open(iName.c_str(), std::ios::trunc);
        else
          outfile.open(iName.c_str(), std::ios::app);
        if (outfile){
        outfile << std::string("# ") + iName<<'\n';
        for (int i =0; i<x.size();i++)
            outfile<< x[i]<<" "<<y[i]<<'\n';
        outfile.flush();
        outfile.close();
        }
        }
        else std::cout<<"Sizes are not equal!";
}


 template<typename X, typename Y,typename Z>
    void writeToFile (std::string iName,const  X& x,const  Y& y,const  Z& z, bool rewrite){
             std::ofstream outfile;
        if (x.size()==y.size()&&x.size()==z.size()){
        if (rewrite)
          outfile.open(iName.c_str(), std::ios::trunc);
        else
          outfile.open(iName.c_str(), std::ios::app);
        if (outfile){
        outfile << std::string("# ") + iName<<'\n';
        for (int i =0; i<x.size();i++)
            outfile<< x[i]<<" "<<y[i]<<" "<<z[i]<<'\n';
        outfile.flush();
        outfile.close();
        }
        }
        else std::cout<<"Sizes are not equal!";
}

#endif
