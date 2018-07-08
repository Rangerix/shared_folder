#include <bits/stdc++.h>
using namespace std;

struct stride {
	string id,pdbseq;
	int category;
};

vector<int> findremove(string s1,int s,int e,int category){
	int i,j,k,l;
    int count=0;
    vector<int> startend;
    vector<vector<char> > alpha;
    for(i=0;i<s1.length();){
        vector<char> temp;
        temp.push_back(s1[i]);
        int l=i+1;
        while(s1[l]=='('){
            temp.push_back(s1[l+1]);
            l+=3;
        }
        i=l;
        alpha.push_back(temp);
    }
    string str="";
	if(category == 1 || category ==3||category==-1){
	    str="";
	    for(i=0;i<alpha.size();i++) 
	    	str+=alpha[i][0];
	    str[s-1]='#';//flag
	    str[e-1]='#';

	    string tempstr="";
	    for(i=0;i<str.length();i++){
	    	if(str[i]=='*') continue;
	    	tempstr+=str[i];
	    }
	    str=tempstr;
	}
	else if(category==2){
		str="";
	    for(i=0;i<alpha.size();i++) 
	    	str+=alpha[i][0];
	    string tempstr="";
	    for(i=0;i<str.length();i++){
	    	if(str[i]=='*') continue;
	    	tempstr+=str[i];
	    }
	    str=tempstr;
	    str[s-1]='#';//flag
	    str[e-1]='#';
	}
    count=0;
    for(i=0;str[i]!='#';i++){
    	count++;
    }
    startend.push_back(count);
    count=0;
    for(i=str.length()-1;str[i]!='#';i--){
    	count++;
    }
    startend.push_back(count);
    return startend;
}

string makepdbfile(string pdbid,int s,int e){
	char chainid=pdbid[5];
	char pdbidch[5];
	int i,j,k,l;
	for(i=0;i<4;i++){
		pdbidch[i]=pdbid[i];
	}
	pdbidch[i++]='\0';
	string filename="../../pdb/pdb";													//pdb file path
	string tempstr(pdbidch);
	filename+=tempstr;
	filename+=".ent";
	transform(filename.begin(),filename.end(),filename.begin(),::tolower);
	char tempfilenamech[]="tempfile.txt";
	string tempfilename(tempfilenamech);
	string str;
	ifstream ifs(filename.c_str());
	if(!ifs){
		printf("problem in opening %s\n",filename.c_str());
		return "";
	}
	transform(filename.begin(),filename.end(),filename.begin(),::tolower);
	ofstream ofs(tempfilenamech);
	while(ifs){
		getline(ifs,str);
		if(!ifs) break;
		if(str.substr(0,6)=="ENDMDL") break;
		ofs<<str<<endl;
	}
	ifs.close();
	ofs.close();
	string command = "grep ^ATOM "+tempfilename+" > mainfile.txt";
    system(command.c_str());

    //extract only the CA LINES
    char s1[]="mainfile.txt";
    char s2[]="mainCAfile.txt";
    ifs.open(s1);
    ofs.open(s2);
    vector<string> alpha;
    while(ifs){
    	getline(ifs,str);
    	if(!ifs) break;
    	if(str[21]==chainid){
    		alpha.push_back(str);
    		ofs<<str<<endl;
    	}
    }
    ifs.close();
    ofs.close();
    str="";
    int count=0;
    string oldnum=alpha[0].substr(22,5),newnum="";
    for(i=1;i<alpha.size()&&count<s;i++){
    	str=alpha[i];
    	newnum=str.substr(22,5);
    	if(newnum!=oldnum){
    		count++;
    	}
    	oldnum=newnum;
    }

    count=0;
    oldnum=alpha[alpha.size()-1].substr(22,5);
    for(j=alpha.size()-2;j>i&&count<e;j--){
    	str=alpha[j];
    	newnum=str.substr(22,5);
    	if(newnum!=oldnum){
    		count++;
    	}
    	oldnum=newnum;
    }

    filename="";
    filename=pdbid;
    filename+=".pdb";
    //cout<<filename<<endl;
    ofs.open(filename.c_str());
    for(k=i;k<=j;k++){
    	ofs<<alpha[k]<<endl;
    }
    return filename;
}

int main(int argc,char* argv[])
{
	if(argc!=3){
		cout<<"exec stridefile blastfilename\n";
		return 1;
	}


	ifstream ifs(argv[1]);
	if(!ifs){
		cout<<"problem "<<argv[1]<<endl;
		return 1;
	}
	vector<struct stride> stridevec;
	struct stride temp;
	int i,j;
	string str;
	while(ifs){
		str="";
		getline(ifs,str);
		if(!ifs) break;
		if(str.length()>8&&str.substr(0,8)=="PDB_ID :"){
			temp.id=str.substr(9,6);
		}
		else if(str.length()>12&&str.substr(0,12)=="PDB sequence"){
			getline(ifs,str);
			if(!ifs) break;
			temp.pdbseq=str;
		}
		else if(str.length()>13&&str.substr(0,13)=="dist sequence"){
			getline(ifs,str);
			if(!ifs) break;
			temp.pdbseq=str;
		}
		else if(str.length()>10&&str.substr(0,5)=="match"){
			if(str[8]=='1')	temp.category=1;							//match : 1
			else if(str[8]=='-'&&str[9]=='1') temp.category=-1;
			else if(str[24]=='1') temp.category=2;						//match_blank : 1
			else if(str[37]=='1') temp.category=3;						
		}
		else if(str==""){
			stridevec.push_back(temp);
			temp.id="";
			temp.pdbseq="";
			temp.category=0;
		}
	}
	ifs.close();

	cout<<"---------------------------- all sequence details stored ---------------------------------\n";

	ifs.open(argv[2]);
	string str1,str2,str3;
	while(ifs){
		str1="";
		getline(ifs,str1);
		if(!ifs) break;
		if(str1[0]=='>'){
			string targetid="";
			for(i=0;str1[i]!=':';i++){
				targetid+=str1[i];
			}
			targetid=targetid.substr(1,6);
			i++;
			string matchid="";
			for(;str1[i]!=':';i++){
				matchid+=str1[i];
			}
			matchid=matchid.substr(0,6);
			int tarstart,tarend,matchstart,matchend;
			i++;
			char tempch[100];
			int k=0;
			for(;str1[i]!=' ';i++)
				tempch[k++]=str1[i];
			tempch[k++]='\0';
			tarstart=atoi(tempch);
			i++;

			k=0;
			for(;str1[i]!=':';i++)
				tempch[k++]=str1[i];
			tempch[k++]='\0';
			tarend=atoi(tempch);
			i++;

			k=0;
			for(;str1[i]!=' ';i++)
				tempch[k++]=str1[i];
			tempch[k++]='\0';
			matchstart=atoi(tempch);
			i++;

			k=0;
			for(;i<str1.length();i++)
				tempch[k++]=str1[i];
			tempch[k++]='\0';
			matchend=atoi(tempch);
			i++;

			//cout<<targetid<<" "<<matchid<<endl;
			for(i=0;i<stridevec.size();i++){
				if(targetid==stridevec[i].id) break;
			}
			//cout<<"here1"<<endl;
			for(j=0;j<stridevec.size();j++){
				if(matchid==stridevec[j].id) break;
			}
			//cout<<"here2"<<endl;
			if(i==stridevec.size()||j==stridevec.size()) continue;
			string tpdbstring=stridevec[i].pdbseq;
			string mpdbstring=stridevec[j].pdbseq;
			vector<int> startendt=findremove(tpdbstring,tarstart,tarend,stridevec[i].category);
			vector<int> startendm=findremove(mpdbstring,matchstart,matchend,stridevec[j].category);

			string tpdbfile=makepdbfile(targetid,startendt[0],startendt[1]);
			//cout<<"here3"<<endl;
			string mpdbfile=makepdbfile(matchid,startendm[0],startendm[1]);

			string foldername=targetid+"_"+matchid;
			cout<<foldername<<endl;
			system(("mkdir "+foldername).c_str());
			system(("mv "+tpdbfile+" "+foldername+"/").c_str());
			system(("mv "+mpdbfile+" "+foldername+"/").c_str());
		}
	}

return 0;
}
