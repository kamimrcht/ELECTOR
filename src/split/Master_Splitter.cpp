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


void updateK(kmer& min, char nuc,int k){
    min<<=2;
    min+=nuc2int(nuc);
    min%=(1<<(2*k));
}



void updateRCK(kmer& min, char nuc,int k){
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
        if(get<0>(next)-get<0>(anchor_start)<1000 and get<0>(next) > get<0>(anchor_start) ){
            if(get<1>(next)-get<1>(anchor_start)<1000 and get<1>(next) > get<1>(anchor_start) ){
                if(get<2>(next)-get<2>(anchor_start)<1000 and  get<2>(next) > get<2>(anchor_start)){
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
    return res/2;
}

string generate_dumb_str(uint n,const string& header,const string& begin,const string& end){
	string res;
	if(not end.empty()){
		res+=header+"\n"+end+"\n";
	}
	for(uint i(0);i<n-1;++i){
		res+=header+"\nN\n";
	}
	if(not begin.empty()){
		res+=header+"\n"+begin+"\n";
	}
	if(begin.empty() and end.empty()){
		res+=header+"\nN\n";
	}
	return res;
}



uint largest_fragment(const string& str){
    uint res(0);
    uint last(0);
    for(uint i(0);i<str.size();++i){
        if(str[i]=='\n'){
            uint len=i-last;
            res=max(res,len);
            last=i;
        }
    }
    return res;
}





void split(const string& ref, const string& S1, const string& S2, string& out_ref, string& out_S1, string& out_S2,const string& header,bool first_call,int k,uint minSize=20){
    unordered_map<kmer,position> kmer_ref,kmer_ref_inS1,kmer_shared;
    kmer seq(str2num(ref.substr(0,k)));
    kmer_ref[seq]=0;
    for(uint j(0);j+k<ref.size();++j){
        updateK(seq,ref[j+k],k);
        if(kmer_ref.count(seq)==0){
            kmer_ref[seq]=j+1;
        }else{
            //repeated kmer in ref
            kmer_ref[seq]=-1;
        }
    }

    seq=str2num(S1.substr(0,k));
    if(kmer_ref.count(seq)==1){
        if(kmer_ref[seq]!=-1){
            kmer_ref_inS1[seq]=0;
        }
    }

    for(uint j(0);j+k<S1.size();++j){
        updateK(seq,S1[j+k],k);
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


    seq=str2num(S2.substr(0,k));
    if(kmer_ref_inS1.count(seq)==1){
        if(kmer_ref_inS1[seq]!=-1){
            kmer_shared[seq]=0;
        }
    }

    for(uint j(0);j+k<S2.size();++j){
        updateK(seq,S2[j+k],k);
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
        updateK(seq,ref[j+k],k);
        if(kmer_shared.count(seq) ){
        if(kmer_shared[seq]!=-1 and j-last_indexed_anchor> minSize){
                anchor_list.push_back(make_tuple(kmer_ref[seq],kmer_ref_inS1[seq],kmer_shared[seq]));
                last_indexed_anchor=j;
            }
        }
    }

    //Anchors list filled Now to find maximal chain
    auto BL(best_chain_from_anchor_list(anchor_list));

    if(BL.size()<1){
		out_ref+=header+"\n"+ref+"\n";
		out_S2+=header+"\n"+S2+"\n";
		out_S1+=header+"\n"+S1+"\n";
		return;
	}
	int i=(0);
    uint pred_ref(0),pred_S1(0),pred_S2(0);
    string start_ref(ref.substr(pred_ref,get<0>(anchor_list[BL[i]])+k));
    string start_S1(S1.substr(pred_S1,get<1>(anchor_list[BL[i]])+k));
    string start_S2(S2.substr(pred_S2,get<2>(anchor_list[BL[i]])+k));
    string out_ref_2,out_S1_2,out_S2_2;
    if(start_S2.size()*2<start_ref.size() and start_ref.size()-start_S2.size()>200 and first_call and true){
		split(start_ref,start_S1,start_ref,out_ref_2,out_S1_2,out_S2_2,header,false,k,1.2*start_S2.size());
		out_ref+=out_ref_2;
		out_S1+=out_S1_2;
		out_S2+=generate_dumb_str(fragment(out_ref_2),header,start_S2,"");
		pred_S1=get<1>(anchor_list[BL[i]])+k;
		pred_ref=get<0>(anchor_list[BL[i]])+k;
		pred_S2=get<2>(anchor_list[BL[i]])+k;
		++i;
	}


    for(;i<(int)BL.size()-1;++i){
        int size_R(get<0>(anchor_list[BL[i]])-pred_ref),size_S1(get<1>(anchor_list[BL[i]])-pred_S1),size_S2(get<2>(anchor_list[BL[i]])-pred_S2);
        if(size_R>minSize and size_S1>minSize and size_S2>minSize and abs(size_S1-size_R)<size_R*0.5 and abs(size_S2-size_R)<size_R*0.5 ){
            out_ref+=header+"\n"+ref.substr(pred_ref,get<0>(anchor_list[BL[i]])-pred_ref+k)+"\n";
            out_S2+=header+"\n"+S2.substr(pred_S2,get<2>(anchor_list[BL[i]])-pred_S2+k)+"\n";
            out_S1+=header+"\n"+S1.substr(pred_S1,get<1>(anchor_list[BL[i]])-pred_S1+k)+"\n";
            pred_S1=get<1>(anchor_list[BL[i]])+k;
            pred_ref=get<0>(anchor_list[BL[i]])+k;
            pred_S2=get<2>(anchor_list[BL[i]])+k;
        }
    }

    string end_ref(ref.substr(pred_ref));
    string end_S1(S1.substr(pred_S1));
    string end_S2(S2.substr(pred_S2));
    if(end_S2.size()*2<end_ref.size() and end_ref.size()-end_S2.size()>200 and first_call and true){
		out_ref_2=out_S1_2=out_S2_2="";
		split(end_ref,end_S1,end_ref,out_ref_2,out_S1_2,out_S2_2,header,false,k,1.2*end_S2.size());
		out_ref+=out_ref_2;
		out_S1+=out_S1_2;
		out_S2+=generate_dumb_str(fragment(out_ref_2),header,"",end_S2);

	}else{
		out_ref+=header+"\n"+ref.substr(pred_ref)+'\n';
		out_S1+=header+"\n"+S1.substr(pred_S1)+'\n';
		out_S2+=header+"\n"+S2.substr(pred_S2)+'\n';
	}

}

void best_split(const string& ref, const string& S1, const string& S2, string& s_ref, string& s_S1, string& s_S2,const string& header){
    int k=15;
    split(ref,S1,S2,s_ref,s_S1,s_S2,header,true,k);
    uint largest_frag(largest_fragment(s_ref));
    while(true){
        k-=2;
        if(k<9){
            return;
        }
        string s_ref_aux,s_S1_aux,s_S2_aux;
        split(ref,S1,S2,s_ref_aux,s_S1_aux,s_S2_aux,header,true,k);
        uint largest_frag_aux(largest_fragment(s_ref_aux));
        if(largest_frag_aux<largest_frag){
            largest_frag=largest_frag_aux;
            s_ref=s_ref_aux;
            s_S1=s_S1_aux;
            s_S2=s_S2_aux;
            s_ref_aux=s_S1_aux=s_S2_aux="";
        }else{
            return;
        }
    }
}


uint64_t count_lines(const string& str){
	uint64_t res(0);
	ifstream inFile(str);
	string lineBuffer;
	while(!inFile.eof())
	{
		//~ cout<<"COUNT"<<endl;
		getline (inFile, lineBuffer);
		getline (inFile, lineBuffer);
		res++;
	}
	//~ cout<<"res:"<<res<<endl;
	return res;
}



int main(int argc, char ** argv){

    string inputRef(argv[1]);
    string inputS1(argv[2]);
    string inputS2(argv[3]);
    string outputRef(argv[4]);
    string outputS1(argv[5]);
    string outputS2(argv[6]);
    int k=(stoi(argv[7]));
    int nb_file=(stoi(argv[8]));
    uint64_t max_nuc_amount=(stoll(argv[9])),nuc_amount(0);
    double SIZE_CORRECTED_READ_THRESHOLD=(stod(argv[10]));
    string outDir(argv[11]);
    int64_t factor((count_lines(inputRef)/(nb_file)));
    //~ cout<<count_lines(inputRef)<<" "<<nb_file<<" "<<factor<<endl;

    factor+=1;
    int small_reads(0);
    int wrong_reads(0);
    string progress_file(outDir + "/progress.txt");
    string ref,S1,S2;
    string href,hS1,hS2,s_ref,s_S1,s_S2,line;
    uint64_t position_ref(0),position_cor(0),position_err(0);
    ifstream inR(inputRef),in1(inputS1),in2(inputS2),progress_in(progress_file);
    if(progress_in.good() and not progress_in.eof()){
        getline(progress_in,line);
        position_ref=stoll(line);
        getline(progress_in,line);
        position_cor=stoll(line);
        getline(progress_in,line);
        position_err=stoll(line);
        inR.seekg (position_ref, inR.beg);
        in1.seekg (position_cor, in1.beg);
        in2.seekg (position_err, in2.beg);
    }

    vector<ofstream> outR(nb_file),out1(nb_file),out2(nb_file);
    for(uint i(0);i<nb_file;++i){
        outR[i].open(outputRef+to_string(i),ofstream::trunc);
        out1[i].open(outputS1+to_string(i),ofstream::trunc);
        out2[i].open(outputS2+to_string(i),ofstream::trunc);
    }
    uint64_t i(0);
    while(not inR.eof() and not in2.eof() and not in1.eof()){
        if(nuc_amount>max_nuc_amount){
            break;
        }
        //~ #pragma omp parallel for ordered schedule(dynamic)
        for(uint ii=(0);ii<1000;++ii){
            if(nuc_amount>max_nuc_amount){
                continue;
            }
            #pragma omp ordered
            {
                getline(inR,href);
                getline(inR,ref);
                getline(in1,hS1);
                getline(in1,S1);
                getline(in2,hS2);
                getline(in2,S2);
            }
            if(ref.size()>2){
				if ((double) S2.size() / ref.size() >= SIZE_CORRECTED_READ_THRESHOLD){
                    best_split(ref,S1,S2,s_ref,s_S1,s_S2,href);
                    if((fragment(s_ref))<=1){
						s_ref=href+"\nAAA\n";
						s_S1=href+"\nAAA\n";
						s_S2=href+"\nAAA\n";
                        #pragma omp atomic
                        wrong_reads++;
                    }else{
                    }
                }else{
					s_ref=href+"\nAAA\n";
					s_S1=href+"\nAAA\n";
					s_S2=href+"\nAAA\n";
                    #pragma omp atomic
                    small_reads++;
                }
				#pragma omp ordered
				{
					//~ cout<<i<<" "<<factor<<" "<<i/factor<<endl;
					outR[i/factor]<<s_ref;
					out1[i/factor]<<s_S1;
					out2[i/factor]<<s_S2;
					nuc_amount+=s_ref.size();
				}
				href=s_ref=s_S1=s_S2=ref=S1=S2="";
				++i;
            }else{
			}

        }
    }
    for(uint i(0);i<nb_file;++i){
        outR[i].close();
        out1[i].close();
        out2[i].close();
    }
    ofstream out_small(outDir + "small_reads.txt");
    ofstream out_wrong(outDir + "wrongly_cor_reads.txt");
    out_small<<small_reads<<endl;
    out_wrong<<wrong_reads<<endl;
    out_small.close();
    out_wrong.close();


    if(inR.eof() or in2.eof() or in1.eof() ){
//	const char* cProgress = (outDir + "/progress.txt").c_str();
//        remove(cProgress);
        return 0;
    }
    ofstream out(progress_file);

    out<<(uint64_t)inR.tellg()<<"\n";
    out<<(uint64_t)in1.tellg()<<"\n";
    out<<(uint64_t)in2.tellg()<<"\n"<<flush;
    out.close();
    return 1;
}
