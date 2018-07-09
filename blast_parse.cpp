#include <bits/stdc++.h>
using namespace std;

string getmatchidf(string str){
	// gnl|num|id => id_num
	string str1="",str2="";
	int i,j,l;
	for(i=0;str[i]!='|';i++) ;
	i++;
	for(;str[i]!='|';i++)
		str1+=str[i];
	i++;
	for(;i<str.length();i++)
		str2+=str[i];
	return (str2+"_"+str1);
}


int main(int argc,char* argv[])
{
	if(argc!=3){
		cout<<"exec blastfile resultfile\n";
		return 1;
	}
	ifstream ifs(argv[1]);
	string str="";
	ofstream ofs;
	ofs.open(argv[2],ios::app);
	if(!ofs){
		cout<<"problem "<<argv[2]<<endl;
		return 1;
	}
	int i;
	while(ifs){
		getline(ifs,str);
		if(!ifs) break;

		if(str.substr(0,7)=="unnamed"){
			stringstream ss(str);
			string targetidf,matchidef,tempstr;
			string target="",match="";
			int tarlen,matlen,tst,ten,mst,men;
			float temp;
			ss>>tempstr;
			ss>>matchidef;
			ss>>tarlen>>matlen>>tst>>ten>>mst>>men;
			ss>>temp;
			ss>>temp;
			ss>>temp;
			ss>>temp;
			ss>>target>>match;
			//check for atleast one blank
			int blank=0;
			for(i=0;i<target.length();i++) 
				if(target[i]=='-') {
					blank=1;
					break;
				}
			for(i=0;i<match.length();i++){
				if(match[i]=='-'){
					blank=1;
					break;
				}
			}
			if(blank==0) continue;
			//coverage cutoff
			if((float)(men-mst)/matlen < 0.5) continue;
			if((float)(ten-tst)/tarlen < 0.5) continue;
			//cout<<target<<endl<<match<<endl;
			//check for identity cutoff
			int total=0,count=0;
			for(i=0;i<target.length();i++){
				if(target[i]=='-'||match[i]=='-') continue;
				total++;
				if(target[i]==match[i]) count++;
			}
			//cout<<count<<" "<<total<<endl;
			if((float)(count)/total<1.0) continue;
			//cout<<target<<endl<<match<<endl;
			targetidf="";
			for(i=0;argv[1][i]!='.';i++) targetidf+=argv[1][i];

			ofs<<">"<<targetidf<<":"<<getmatchidf(matchidef)<<":"<<tst<<" "<<ten<<":"<<mst<<" "<<men<<endl;
			ofs<<target<<endl;
			ofs<<match<<endl;

		}
	}
	ifs.close();
	ofs.close();
	

return 0;
}