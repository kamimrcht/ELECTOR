#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>



using namespace std;


int main(int argc, char ** argv){
	if(argc<3){
		cout<<"[msa file in]  [msa file out]"<<endl;
		exit(0);
	}

	string input(argv[1]),output(argv[2]),ref, useless,acc1,acc2,acc3,header,cor,err;
	srand (time(NULL));
	ifstream in(input);
	ofstream out(output,ofstream::app);
	if(not in){
		cout<<"Problem opening file"<<endl;
		return 0;
	}

	getline(in,useless);
	getline(in,ref);
	getline(in,useless);
	getline(in,err);
	getline(in,useless);
	getline(in,cor);
	acc1+=ref;
	acc2+=err;
	acc3+=cor;
	header=useless;

	while(not in.eof()){
		getline(in,useless);
		getline(in,ref);
		getline(in,useless);
		getline(in,err);
		getline(in,useless);
		getline(in,cor);
		if(header!=useless){
			if(acc1.size()>1){
				out<<header.substr(0, header.size()-11)<<" \n"<<acc1<<'\n';
				out<<header.substr(0, header.size()-11)<<" \n"<<acc2<<'\n';
				out<<header.substr(0, header.size()-11)<<" \n"<<acc3<<'\n';
				header=useless;
			}
			acc1=ref;
			acc2=err;
			acc3=cor;
		}else{
			acc1+=ref;
			acc2+=err;
			acc3+=cor;
		}
	}
	//~ if(acc1.size()>1){
		//~ out<<header<<'\n'<<acc1<<'\n';
		//~ out<<header<<'\n'<<acc2<<'\n';
		//~ out<<header<<'\n'<<acc3<<'\n';
		//~ header=useless;
		//~ acc1=acc2=acc3="";
	//~ }
	out.close();
}
