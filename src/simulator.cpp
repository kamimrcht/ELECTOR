#include <fstream>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>



using namespace std;


static unsigned int seed;


uint32_t xs(uint32_t& y){
	y^=(y<<13); y^=(y>>17);y=(y^=(y<<15)); return y;
}





char randNucle(char c='N'){
	//~ switch (rand()%4){
	switch (xs(seed)%4){
		case 0:
			if(c!='A'){
				return 'A';
			}
			return randNucle(c);
		case 1:
			if(c!='C'){
				return 'C';
			}
			return randNucle(c);
		case 2:
			if(c!='G'){
				return 'G';
			}
			return randNucle(c);
		case 3:
			if(c!='T'){
				return 'T';
			}
			return randNucle(c);
	}
	return randNucle(c);
}


void insertion(double rate, string& result){
	uint dice(rand() % 100);
	if(dice < rate){
		char newNucleotide(randNucle());
		result.push_back(newNucleotide);
		insertion(rate, result);
	}
}


string mutateSequence(const string& referenceSequence,uint mutRate, vector <double> ratioMutation={0.06,0.73,0.21}){
	string result;
	result.reserve(5 * referenceSequence.size());
	for(uint i(0); i < referenceSequence.size(); ++i){
		double substitutionRate(mutRate * ratioMutation[0]);
		double insertionRate(mutRate * ratioMutation[1]);
		double deletionRate(mutRate * ratioMutation[2]);
		uint dice(rand() % 100);


		if (dice <substitutionRate ){
			//SUBSTITUTION
			char newNucleotide(randNucle());
			while(newNucleotide == referenceSequence[i]){
				newNucleotide = randNucle();
			}
			result.push_back(newNucleotide);
			continue;
		} else if(dice < deletionRate+substitutionRate){
			//DELETION
			uint dice2(rand() % 100);
			while (dice2 < deletionRate+substitutionRate){ // deletions larger than 1
				++i;
				dice2 = rand() % 100;
			}
			continue;
		} else if (dice < deletionRate + substitutionRate + insertionRate){
			//INSERTION
			char newNucleotide(randNucle());
			result.push_back(referenceSequence[i]);
			result.push_back(newNucleotide);
			insertion(deletionRate + substitutionRate + insertionRate, result); // larger than 1 insertions

			continue;
		} else {
			result.push_back(referenceSequence[i]);
		}

	}
	return result;
}







int main(int argc, char ** argv){
	if(argc<5){
		cout<<"[Genome reference file] [read length] [coverage] [error rate] [prefix] [LR]"<<endl;
		exit(0);
	}
	bool long_reads(false);
	//~ if(argc==6){
		long_reads=true;
	//~ }
	string input(argv[1]);
	double coverage(stof(argv[3]));
	float length(stof(argv[2]));
	srand (time(NULL));
	ifstream in(input);
	uint errorRate((stof(argv[4]))*10000);
	string prefix(argv[5]);
	string useless, ref,read,pread;
	uint i(0);
	ofstream perfect("p."+prefix+".fa"),out(prefix+".fa");
	while(not in.eof()){
		getline(in,useless);
		getline(in,ref);
		if(not ref.empty() and not useless.empty()){
			uint64_t nucProduced(0);
			while(nucProduced<(uint64_t)(coverage*ref.size())){
				if(i%100==0){
					seed=(rand());
				}
				//produce a read
				uint64_t position=xs(seed)%ref.size();
				//~ uint position=rand()%ref.size();
				if(position+length<=ref.size()){
					bool valid(true);
					uint error(0);
					pread=ref.substr(position,length);
					read=pread;
					if(long_reads){
						read=mutateSequence(read,errorRate/100);
					}else{
						for(uint i(0);i<read.size();++i){
							if(read[i]=='N' or read[i]=='n'){valid=false;break;}
							if(xs(seed)%10000<=errorRate){
							//~ if(rand()%errorRate==0){
								read[i]=randNucle(read[i]);
								++error;
							}
						}
					}
					if(valid){
						perfect<<">"+to_string(i)<<"\n";
						perfect<<pread<<"\n";
						out<<">"+to_string(i)<<"\n";
						out<<read<<"\n";
						nucProduced+=read.size();
						++i;
					}
				}
			}

		}
	}
}
