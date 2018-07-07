#include <bits/stdc++.h>
using namespace std;

int main(int argc,char* argv[])
{
	if(argc!=2){
		cout<<"./a.out sequence-file\n";
		return 1;
	}
	ifstream ifs(argv[1]);
	string str1="",str2="";
	int k,i,j,l;
	string filename="";
	while(ifs){
		getline(ifs,str1);
		if(!ifs) break;
		long long int num=-1;
		string chnum="";
		if(str1[0]=='>'){
			cout<<filename<<endl;
			ofstream ofs((filename+".seq").c_str());
			ofs<<str2;
			ofs.close();
			string command="/home/kkg/Documents/databases/ncbi-blast-2.2.30+/bin/blastp -query ";
			command+=filename; 
			command+=".seq -db /home/kkg/Documents/databases/ncbi-blast-2.2.30+/bin/allseq -outfmt \"7 qseqid sseqid qlen slen qstart qend sstart send evalue pident gaps score qseq sseq\" -out ";
			command+=filename;
			command+=".blast -evalue 0.00001";
			//cout<<command<<endl;
			system(command.c_str());

			k=0;
			chnum="";
			filename="";
			for(i=5;str1[i]!='|';i++){
				chnum+=str1[i];
			}
			i++;
			for(;i<str1.length();i++){
				filename+=str1[i];
			}
			filename+="_";
			filename+=chnum;
			//cout<<"filename : "<<filename<<endl;
			str2="";
		}
		else str2+=str1;
		if(!ifs) break;
		
	}

return 0;
}
