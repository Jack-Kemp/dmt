#ifndef INPUT_OUTPUT_UTILITIES_H
#define INPUT_OUTPUT_UTILITIES_H


#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<map>
#include<functional>
#include<vector>
#include<filesystem>

template<class A>
std::string to_string_sstream(const A& num){	
  std::ostringstream streamObj;
  streamObj << num;
  return streamObj.str();
}


bool findCheckpointFile(std::string checkpointName,
			std::string dir,
			std::string ext,
			double& tStart,
			std::string& checkDMTName)
{
  std::string checkPath;
  bool found = false;
  uint clen = checkpointName.size();
  uint elen = ext.size();
  for (const auto & entry : std::filesystem::directory_iterator(dir))
	{
	  auto path = entry.path().filename().string();
	  if (path.substr(0, clen) == checkpointName
	      and path.substr(path.size()-elen, 4) == ext)
	    {
	      if (path.substr(clen, 3) != "_t_")
		{
		std::cout<<"Checkpoint name match, but no time '_t_' suffix; ignoring."
			 << std::endl;
		}
	      else
		{
		  found = true;
		  auto tCheck = std::stod(path.substr(clen+3, path.size()-clen-elen-3));
		  if (tCheck > tStart)
		    {
		      tStart = tCheck;
		      checkDMTName = dir + '/' + path;
		    }
		}
	    }
	}
  return found;
}

//Writes square data tabulated space delimited to a file, with column
//headers as a comment on top and run information and argument file
//as comments at the end.
template<class X, class Y, class Z>
void writeDataToFile (std::string outName, const X& data, const Y& data2D, const Z& runInfo, std::string argFileName){
  std::ifstream argFile;
  argFile.open(argFileName.c_str());
  std::ofstream outfile;
  outfile.open(outName.c_str(), std::ios::trunc);
  std::vector<std::string> ignoreCols;

  //Check data is actually square, ignore cols with nrows < max(nrows)
  uint size = 0;
  for (auto & [key, value] : data2D)
    size = value.size() > size ? value.size() : size;
  for (auto & [key, value] : data)
    size = value.size() > size ? value.size() : size;
  for (auto & [key, value] : data2D)
    if(size != value.size()){
      std::cout<<key<<" has only "<< value.size() << "rows, not writing.";
      ignoreCols.push_back(key);
    }
  for (auto & [key, value] : data)
    if(size != value.size()){
      std::cout<<key<<" has only "<< value.size() << " rows, not writing.";
      ignoreCols.push_back(key);
    }
  
  if(outfile and argFile){
    //Write the column labels
    outfile << '#';
    for (auto & [key, value] : data2D)
      if(std::find(ignoreCols.begin(), ignoreCols.end(), key) == ignoreCols.end())
	for(uint i = 0; i < value[0].size(); i++)
	  outfile << key +  '_' + to_string_sstream(i) + ' ';
    for (auto & [key, value] : data)
      if(std::find(ignoreCols.begin(), ignoreCols.end(), key) == ignoreCols.end())
	outfile << key + ' ';
    outfile << '\n';

    //Write the data in rows
    for (uint i = 0; i<size; i++){
      for (auto & [key, value] : data2D)
	if(std::find(ignoreCols.begin(), ignoreCols.end(), key) == ignoreCols.end())
	  for(uint j = 0; j < value[0].size(); j++)
	    outfile << value[i][j] << " ";
      for(auto & [key, value] : data)
	if(std::find(ignoreCols.begin(), ignoreCols.end(), key) == ignoreCols.end())
	  outfile << value[i] << ' ';
      outfile << '\n';
    }
    
    outfile <<"\n###############################################\n";
    outfile <<"#                Run Information\n";
    outfile <<"###############################################\n";
    for(auto & [key, value] : runInfo)
      outfile <<"#:" << key <<" = " << value <<"\n";
    outfile <<"\n###############################################\n";
    outfile <<"#                #Input Parameters\n";
    outfile <<"###############################################\n";
    std::string str; 
    while (std::getline(argFile, str))
      {
	outfile << "#~" + str + '\n';
      }
    outfile.flush();
    outfile.close();
  }
  else std::cout<<"Could not open file!"<<std::endl;
}




template <class A, class B, class C>
std::string construct_label(const std::string& base,
                            const std::string& vary,
                            A reals, B bools, C ints){
  std::string label = base;
  if (bools["Vectorize"]) label +=  "_vec";
  if (bools["HermitianBasis"]) label +=  "_herm";
  if (bools["ConserveQNs"]) label += "_cons";
      
  for (auto const& x : reals)
    {
      label += "_" + x.first + "_";
      label += (x.first == vary ? "vary" : to_string_sstream(x.second));
    }
  for (auto const& x : ints)
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
    for  (uint i =0; i<x.size()&& file >> x[i]; i++);
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
    for (uint i =0; i<x.size();i++)
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
      for (uint i =0; i<x.size();i++)
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
      for (uint i =0; i<x.size();i++)
	outfile<< x[i]<<" "<<y[i]<<" "<<z[i]<<'\n';
      outfile.flush();
      outfile.close();
    }
  }
  else std::cout<<"Sizes are not equal!";
}

template<typename X>
void write2DToFile (std::string iName,const  X& x, bool rewrite = true){
  if(x.size() > 0){
    if(x[0].size() > 0){
      std::ofstream outfile;
      if (rewrite)
	outfile.open(iName.c_str(), std::ios::trunc);
      else
	outfile.open(iName.c_str(), std::ios::app);
      if (outfile){
	outfile << std::string("# ") + iName<<'\n';
	for (uint i =0; i<x.size();i++)
	  {
	    for(uint j=0; j<x[0].size();j++)
	      outfile<< x[i][j]<<' ';
	    outfile<<'\n';
	  }
	outfile.flush();
	outfile.close();
      }
    }
  }
}

constexpr unsigned int hash(const char *s, int off = 0) {                        
    return !s[off] ? 5381 : (hash(s, off+1)*33) ^ s[off];                           
}                                                                                

constexpr inline unsigned int operator "" _(char const * p, size_t) { return hash(p); }


#endif
