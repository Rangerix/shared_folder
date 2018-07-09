#include <bits/stdc++.h>
using namespace  std;

int main(int argc,char* argv[]){
	if(argc!=3){
		cout<<"exec blastfilelist resultfile\n";
		return 1;
	}
	ifstream ifs(argv[1]);
	string str;
	string filename(argv[2]);
	while(ifs){
		getline(ifs,str);
		if(!ifs) break;
		string command="./exec ";
		command+=str;
		command+=" ";
		command+=filename;
		cout<<command<<endl;
		system(command.c_str());
	}
	ifs.close();


return 0;
}
