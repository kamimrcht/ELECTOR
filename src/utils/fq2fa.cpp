#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>



using namespace std;



int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[Fasta file]"<<endl;
		exit(0);
	}
	string input(argv[1]);
	srand (time(NULL));
	string ref, useless,header;
	ifstream in(input);
	vector<uint> lengths;
	while(not in.eof()){
		getline(in,header);
		getline(in,ref);
		getline(in,useless);
		getline(in,useless);
		if(ref.size()>1){
			cout<<">"+header.substr(1)<<"\n";
			cout<<ref<<"\n";
		}
		ref="";
	}
}
