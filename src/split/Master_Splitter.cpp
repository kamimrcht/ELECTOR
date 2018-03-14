#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <tuple>
#include <algorithm>



using namespace std;



typedef uint32_t kmer;
typedef int32_t position;
typedef tuple<position,position,position> anchor;



int k;



kmer str2num(const string& str){
	kmer res(0);
	for(uint64_t i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}


kmer nuc2int(char c){
	switch(c){
		//~ case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	return 0;
}



kmer nuc2intrc(char c){
	switch(c){
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		//~ case 'T': return 0;
	}
	return 0;
}


void updateK(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=(1<<(2*k));
}



void updateRCK(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*(k-1)));
}



pair<int,int> best_chain_from_anchor(unordered_map<uint,pair<int,int>>& best_chain_computed,const vector<anchor>& anchor_list,uint anchor_indice){
	if(best_chain_computed.count(anchor_indice)==1){
		return best_chain_computed[anchor_indice];
	}

	int max_chain(-1),next_anchor(-1);
	anchor anchor_start(anchor_list[anchor_indice]);
	for(uint i(anchor_indice+1);i<anchor_list.size();++i){
		anchor next(anchor_list[i]);

		if(get<0>(next)-get<0>(anchor_start)<10000 and get<0>(next) > get<0>(anchor_start) ){
			if(get<1>(next)-get<1>(anchor_start)<10000 and get<1>(next) > get<1>(anchor_start) ){
				if(get<2>(next)-get<2>(anchor_start)<10000 and  get<2>(next) > get<2>(anchor_start)){
					auto p=best_chain_from_anchor(best_chain_computed,anchor_list,i);
					if(p.first>max_chain){
						max_chain=p.first;
						next_anchor=i;
					}
				}
			}
		}else{
			//TOO FAR NOW
			break;
		}
	}
	best_chain_computed[anchor_indice]={1+max_chain,next_anchor};
	return {1+max_chain,next_anchor};
}



vector<int> best_chain_from_anchor_list(const vector<anchor>& anchor_list){
	unordered_map<uint,pair<int,int>> best_chain_computed;
	vector<int> res;
	int max_chain(-1),next_anchor(-1);
	for(uint i(0);i<anchor_list.size();++i){
		auto p=best_chain_from_anchor(best_chain_computed,anchor_list,i);
		if(p.first>max_chain){
			max_chain=p.first;
			next_anchor=i;
		}
	}

	while(next_anchor!=-1){
		res.push_back(next_anchor);
		next_anchor=best_chain_computed[next_anchor].second;
	}
	return res;
}


uint fragment(const string& str){
	uint res(0);
	for(uint i(0);i<str.size();++i){
		if(str[i]=='\n'){
			++res;
		}
	}
	return res;
}



void split(const string& ref, const string& S1, const string& S2, string& out_ref, string& out_S1, string& out_S2,const string& header){
	//~ cout<<"GO SPLIT"<<endl;
	unordered_map<kmer,position> kmer_ref,kmer_ref_inS1,kmer_shared;
	kmer seq(str2num(ref.substr(0,k)));
	kmer_ref[seq]=0;
	for(uint j(0);j+k<ref.size();++j){
		updateK(seq,ref[j+k]);
		if(kmer_ref.count(seq)==0){
			kmer_ref[seq]=j+1;
		}else{
			//repeated kmer in ref
			kmer_ref[seq]=-1;
		}
	}

	seq=str2num(S1.substr(0,k));
	if(kmer_ref.count(seq)==1){
		if(kmer_ref[seq]=-1){
			kmer_ref_inS1[seq]=0;
		}
	}

	for(uint j(0);j+k<S1.size();++j){
		updateK(seq,S1[j+k]);
		if(kmer_ref.count(seq)==1){
			if(kmer_ref[seq]!=-1){
				if(kmer_ref_inS1.count(seq)==0){
					kmer_ref_inS1[seq]=j+1;
				}else{
					//repeated kmer
					kmer_ref_inS1[seq]=-1;
				}
			}
		}
	}

	//~ cout<<kmer_ref_inS1.size()<<endl;

	seq=str2num(S2.substr(0,k));
	if(kmer_ref_inS1.count(seq)==1){
		if(kmer_ref_inS1[seq]!=-1){
			kmer_shared[seq]=0;
		}
	}

	for(uint j(0);j+k<S2.size();++j){
		updateK(seq,S2[j+k]);
		if(kmer_ref_inS1.count(seq)==1){
			if(kmer_ref_inS1[seq]!=-1){
				if(kmer_shared.count(seq)==0){
					kmer_shared[seq]=j+1;
				}else{
					//repeated kmer
					kmer_shared[seq]=-1;
				}
			}
		}
	}
	//~ cout<<S2.size()<<endl;
	//~ cout<<"Kmer shared"<<endl;
	//~ if(kmer_shared.size()==0){
		//~ cout<<":("<<flush;
		//~ cout<<ref.size()<<endl;
		//~ cout<<S1.size()<<endl;
		//~ cout<<S2.size()<<endl;
	//~ }
	//~ cout<<kmer_shared.size()<<endl;

	//NOW KMER_SHARED CONTAIN KMER IN COMMUM IN THE THREE READ AND DUPLICATED IN NON OF THE THREE READS

	vector<anchor> anchor_list;
	seq=(str2num(ref.substr(0,k)));
	if(kmer_shared.count(seq)){
		if(kmer_shared[seq]!=-1){
			anchor_list.push_back(make_tuple(kmer_ref[seq],kmer_ref_inS1[seq],kmer_shared[seq]));
		}
	}

	uint last_indexed_anchor(0);
	for(uint j(0);j+k<ref.size();++j){
		updateK(seq,ref[j+k]);
		if(kmer_shared.count(seq) ){
		if(kmer_shared[seq]!=-1 and j-last_indexed_anchor> 10){
				anchor_list.push_back(make_tuple(kmer_ref[seq],kmer_ref_inS1[seq],kmer_shared[seq]));
				last_indexed_anchor=j;
			}
		}
	}

	//Anchors list filled Now to find maximal chain
	auto BL(best_chain_from_anchor_list(anchor_list));

	//~ cout<<"best chain size"<<endl;
	//~ cout<<BL.size()<<endl;
	uint pred_ref(0),pred_S1(0),pred_S2(0);
	for(int i(0);i<(int)BL.size()-1;++i){
		//~ cout<<"NO"<<endl;
		if(get<0>(anchor_list[BL[i]])-pred_ref>15 and get<2>(anchor_list[BL[i]])-pred_S2>15 and get<1>(anchor_list[BL[i]])-pred_S1>15 ){
			out_ref+=header+"\n"+ref.substr(pred_ref,get<0>(anchor_list[BL[i]])-pred_ref+k)+"\n";
			out_S2+=header+"\n"+S2.substr(pred_S2,get<2>(anchor_list[BL[i]])-pred_S2+k)+"\n";
			out_S1+=header+"\n"+S1.substr(pred_S1,get<1>(anchor_list[BL[i]])-pred_S1+k)+"\n";
			pred_S1=get<1>(anchor_list[BL[i]])+k;
			pred_ref=get<0>(anchor_list[BL[i]])+k;
			pred_S2=get<2>(anchor_list[BL[i]])+k;
		}
	}
	//~ cout<<"??"<<endl;
	out_ref+=header+"\n"+ref.substr(pred_ref)+'\n';
	out_S1+=header+"\n"+S1.substr(pred_S1)+'\n';
	out_S2+=header+"\n"+S2.substr(pred_S2)+'\n';
	//~ cout<<"?"<<endl;
}

void best_split(const string& ref, const string& S1, const string& S2, string& s_ref, string& s_S1, string& s_S2,const string& header){
	k=15;
	split(ref,S1,S2,s_ref,s_S1,s_S2,header);
	uint nb_frag(fragment(s_ref));
	while(true){
		k-=2;
		if(k<7){
			return;
		}
		string s_ref_aux,s_S1_aux,s_S2_aux;
		split(ref,S1,S2,s_ref_aux,s_S1_aux,s_S2_aux,header);
		uint nb_frag_aux(fragment(s_ref_aux));
		if(nb_frag_aux>=nb_frag){
			nb_frag=nb_frag_aux;
			s_ref=s_ref_aux;
			s_S1=s_S1_aux;
			s_S2=s_S2_aux;
			s_ref_aux=s_S1_aux=s_S2_aux="";
		}else{
			return;
		}
	}
}



int main(int argc, char ** argv){

	string inputRef(argv[1]);
	string inputS1(argv[2]);
	string inputS2(argv[3]);
	string outputRef(argv[4]);
	string outputS1(argv[5]);
	string outputS2(argv[6]);
	k=(stoi(argv[7]));
	int nb_file=(stoi(argv[8]));
	int max_nuc_amount=(stoi(argv[9])),nuc_amount(0);
	int SIZE_CORRECTED_READ_THRESHOLD=(stoi(argv[10]));
	int skipped_reads(0);

	string progress_file("progress.txt");
	string ref,S1,S2;
	string href,hS1,hS2,s_ref,s_S1,s_S2,line;
	uint32_t position_ref(0),position_cor(0),position_err(0);
	//~ cout<<1<<endl;
	ifstream inR(inputRef),in1(inputS1),in2(inputS2),progress_in(progress_file);
	//~ cout<<"Teenage mutant NINJA TURTLE"<<endl;
	if(progress_in.good() and not progress_in.eof()){
		getline(progress_in,line);
		position_ref=stoi(line);
		getline(progress_in,line);
		position_cor=stoi(line);
		getline(progress_in,line);
		position_err=stoi(line);
		inR.seekg (position_ref, inR.beg);
		in1.seekg (position_cor, in1.beg);
		in2.seekg (position_err, in2.beg);
		//~ cout<<position_ref<<endl;
	}
	//~ cout<<2<<endl;

	vector<ofstream> outR(nb_file),out1(nb_file),out2(nb_file);
	for(uint i(0);i<nb_file;++i){
		outR[i].open(outputRef+to_string(i),ofstream::trunc);
		out1[i].open(outputS1+to_string(i),ofstream::trunc);
		out2[i].open(outputS2+to_string(i),ofstream::trunc);
	}
	//~ cout<<"WTF11"<<endl;
	uint i(0);
	while(not inR.eof() and not in2.eof() and not in1.eof()){
		//~ cout<<4<<endl;
		getline(inR,href);

		getline(inR,ref);
		//~ cout<<ref<<endl;
		getline(in1,hS1);
		getline(in1,S1);
		getline(in2,hS2);
		getline(in2,S2);
		if(ref.size()>2){
			//~ cout<<4<<endl;
			//~ cout<<(double)S2.size()/ref.size()<<endl;
			//~ cout<<SIZE_CORRECTED_READ_THRESHOLD<<endl;
			if((double)ref.size()/S2.size()<=SIZE_CORRECTED_READ_THRESHOLD){
				//~ cout<<5<<endl;
				best_split(ref,S1,S2,s_ref,s_S1,s_S2,href);
				//~ cout<<k<<endl;
				//~ cout<<6<<endl;
				outR[i%nb_file]<<s_ref;
				out1[i%nb_file]<<s_S1;
				out2[i%nb_file]<<s_S2;
				nuc_amount+=s_ref.size();
				if(nuc_amount>max_nuc_amount){
					//~ cout<<7<<endl;
					break;
				}
			}else{
				skipped_reads++;
			}
		}
		s_ref=s_S1=s_S2=ref=S1=S2="";
		++i;
	}
	//~ cout<<"WTF1"<<endl;
	for(uint i(0);i<nb_file;++i){
		outR[i].close();
		out1[i].close();
		out2[i].close();
	}

	if(inR.eof() or in2.eof() or in1.eof() ){
		remove("progress.txt");
		//~ cout<<"I ENDED"<<endl;
		return (1+skipped_reads);
	}
	ofstream out(progress_file);

	//~ cout<<inR.tellg()<<endl;
	out<<inR.tellg()<<"\n";
	out<<in1.tellg()<<"\n";
	out<<in2.tellg()<<"\n"<<flush;
	out.close();
	//~ cout<<"SO CLOSE"<<endl;
	//~ cout<<1<<flush;
	return -(1+skipped_reads);
}
