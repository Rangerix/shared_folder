#include <bits/stdc++.h>
using namespace std;

char shorten(string str){
    map<string,char> m;
    m["CYS"]= 'C'; m["ASP"]= 'D'; m["SER"]= 'S'; m["GLN"]= 'Q'; m["LYS"]= 'K';
     m["ILE"]= 'I'; m["PRO"]= 'P'; m["THR"]= 'T'; m["PHE"]= 'F'; m["ASN"]= 'N'; 
     m["GLY"]= 'G'; m["HIS"]= 'H'; m["LEU"]= 'L'; m["ARG"]= 'R'; m["TRP"]= 'W'; 
     m["ALA"]= 'A'; m["VAL"]='V';m["GLU"]= 'E'; m["TYR"]= 'Y'; m["MET"]= 'M';

    if(m.find(str)==m.end())
        return 'X';
    else return m[str];
}

char structurerepresent(char ch){
	if(ch=='H'||ch=='G'||ch=='I') return 'H';
	if(ch=='E') return 'S';
	if(ch=='b'||ch=='B'||ch=='T'||ch=='C') return 'L';
}

char arearepresent(char res,string val){
	float given=0;
	if(res == 'A')  { given= 106; }
    else if(res == 'B')  { given= 160; }
    else if(res == 'C')  { given= 135; }
    else if(res == 'D')  { given= 163; }
    else if(res == 'E')  { given= 194; }
    else if(res == 'F')  { given= 197; }
    else if(res == 'G')  { given= 84; }
    else if(res == 'H')  { given= 184; }
    else if(res == 'I')  { given= 169; }
    else if(res == 'K')  { given= 205; }
    else if(res == 'L')  { given= 164; }
    else if(res == 'M')  { given= 188; }
    else if(res == 'N')  { given= 157; }
    else if(res == 'P')  { given= 136; }
    else if(res == 'Q')  { given= 198; }
    else if(res == 'R')  { given= 248; }
    else if(res == 'S')  { given= 130; }
    else if(res == 'T')  { given= 142; }
    else if(res == 'V')  { given= 142; }
    else if(res == 'W')  { given= 227; }
    else if(res == 'X')  { given= 180; }
    else if(res == 'Y')  { given= 222; }
    else if(res == 'Z')  { given= 196; }

    float value=atof(val.c_str());
    float ratio=value/given;
    if(ratio<=0.09) return 'B';
    else if(ratio>0.09 && ratio < 0.64) return 'I';
    else if(ratio>=0.64) return 'E';
}

struct stride{
    char res;
    char secstr;
    char area;
    void show(void){
        cout<<res<<" "<<secstr<<" "<<area<<" ";
    }
};

string runstride(string pdbchain){
	int i,j,k,l;
	char pdbidch[5];
	char chainid=pdbchain[5];
	for(i=0;i<4;i++)
		pdbidch[i]=pdbchain[i];
	pdbidch[i]='\0';
	string pdbid(pdbidch);
	transform(pdbid.begin(),pdbid.end(),pdbid.begin(),::tolower);
	string filename="../pdb/pdb";                                                   //specify path for pdb files
	filename+=pdbid;
	filename+=".ent";

	char tempfilenamech[]="tempfile.txt";
	string tempfilename(tempfilenamech);
	string str;
	ifstream ifs(filename.c_str());
	if(!ifs){
		cout<<"problem in opening \n"<<filename<<endl;
		return "";
	}
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
    char s2[]="mainCAfile.pdb";
    ifs.open(s1);
    ofs.open(s2);
    while(ifs){
    	getline(ifs,str);
    	if(!ifs) break;
    	if(/*str[13]=='C'&&str[14]=='A'&&*/str[21]==chainid)
    		ofs<<str<<endl;
    }
    ifs.close();
    ofs.close();
    string strdfilename=pdbchain+".strd";
    command="/home/kkg/Documents/databases/stride/stride mainCAfile.pdb > "+strdfilename; //specify path for the stride file inside stride folder
    cout<<command<<endl;
    system(command.c_str());

    ////////////////////////////////
    return strdfilename;
    
}

void getsequences(string strdfilename,string pdbchain,char* resultfile){
    //cout<<"getsequences "<<strdfilename<<endl;
    ifstream ifs;
    ifs.open(strdfilename.c_str());
    if(!ifs){
        cout<<"stride file not generated\n";
        return;
    }
    string str="";
    while(ifs){
        getline(ifs,str);
        if(!ifs) break;
        if(str.substr(43,5)=="-Phi-"&&str.substr(53,5)=="-Psi-")
            break;
    }

    int pdbidx=0;
    
    vector<struct stride> stridevec;
    while(ifs){
        getline(ifs,str);
        if(!ifs) break;
        struct stride temp;
        temp.res=shorten(str.substr(5,3));
        temp.secstr=structurerepresent(str[24]);
        temp.area=arearepresent(temp.res,str.substr(61,8));
        stridevec.push_back(temp);
    }
    ifs.close();

    int i,j,k,l;
    vector<vector<char> > alpha;
    for(i=0;i<pdbchain.length();){
        vector<char> temp;
        temp.push_back(pdbchain[i]);
        int l=i+1;
        while(pdbchain[l]=='('){
            temp.push_back(pdbchain[l+1]);
            l+=3;
        }
        i=l;
        alpha.push_back(temp);
    }

    /*int index=0;
    for(index=0;index<alpha.size();index++){
        stridevec[index].show();
        for(i=0;i<alpha[index].size();i++)
            cout<<alpha[index][i]<<" ";
        cout<<endl;
    }*/
    string structureseq="";
    string areaseq="";

    int alphaindex=0;
    int strideidx=0;
    for(alphaindex=0;alphaindex<alpha.size();alphaindex++){
        if(alpha[alphaindex][0]=='*'){
            structureseq+='*';
            areaseq+='*';
            continue;
        }
        int flag=0;
        for(i=0;i<alpha[alphaindex].size();i++){
            if(alpha[alphaindex][i]==stridevec[strideidx].res) {
                structureseq+=stridevec[strideidx].secstr;
                areaseq+=stridevec[strideidx].area;
                strideidx++;
                flag=1;
                break;
            }
        }
        if(flag==0) {
            structureseq+=' ';
            areaseq+=' ';
        }
    }
    ofstream ofs;
    ofs.open(resultfile,ios::app);
    ofs<<"secondary structure : \n";
    ofs<<structureseq<<endl;
    ofs<<"solvent accessibility : \n";
    ofs<<areaseq<<endl;
    ofs.close();
}

int main(int argc,char* argv[])
{
	if(argc!=3){
		cout<<"exec alldetailsfile resultfile\n";
		return 1;
	}
	ifstream ifs(argv[1]);
	if(!ifs) {
		cout<<"problem "<<argv[1]<<endl;
		return 1;
	}
    string pdbid;
    string str,strdfilename;
    char notrun[]="stridenotrun.txt";                                   //specify path 
    ofstream ofs,rfs;
    ofs.open(argv[2],ios::app);
	while(ifs){
        str="";
		getline(ifs,str);
		if(!ifs) break;
        if(str!=""&&str.substr(0,6)=="PDB_ID"){
            pdbid=str.substr(9,6);
            strdfilename=runstride(pdbid);
            ifstream strdfs(strdfilename.c_str());
            
            strdfs.close();
            ofs<<str<<endl;
        }
        else if(str.substr(0,14)=="PDB sequence :"){
            ofs<<str<<endl;
            getline(ifs,str);
            ofs<<str<<endl;
            ofs.close();
            getsequences(strdfilename,str,argv[2]);
            ofs.open(argv[2],ios::app);
        }
        else ofs<<str<<endl;
	}
	ofs.close();
return 0;
}