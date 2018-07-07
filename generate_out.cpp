#include <bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[])
{
	if(argc!=3){
		cout<<"./a.out files.txt outputfilename\n";
		return 1;
	}
	ifstream ifs(argv[1]);
	ofstream ofs(argv[2]);
	if(!ifs){
		cout<<"problem in opening"<<argv[1]<<endl;
		return 1;
	}
	string str1="",str2="";
	string pdbid="";
	int count,i,j;
	int flag=0;
	while(ifs){
		getline(ifs,str1);
		if(!ifs) break;
		if(str1.substr(0,6)=="PDB_ID"){
			pdbid=str1.substr(9,6);
			//number=str1.substr(49,52);
			ofs<<">gnl|"<<str1.substr(48)<<"|"<<pdbid<<endl;
		}
		if(str1.substr(0,5)=="fasta"){
			getline(ifs,str2);
			if(!ifs) break;
			j=1;
			for(i=0;i<str2.length();i++){
				ofs<<str2[i];
				flag=0;
				if(j%80==0){
					ofs<<endl;
					flag=1;	
				} 
				j++;
			}
			if(flag==0) ofs<<endl;
		}
		str1="";
		str2="";
		pdbid="";
	}
	ifs.close();
	ofs.close();

return 0;
}