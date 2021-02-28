#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "graph.hpp"
#include "GRN.hpp"
using namespace std;

int main(){
  GRN g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 1);
  bool connected = g.is_DAG();
  // bool connected = true;
  if(connected){
    cout << "Graph is a DAG" << endl;
  }else{
    cout << "Graph is not a DAG" << endl;
  }
}
