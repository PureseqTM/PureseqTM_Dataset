#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector> 
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cmath>
using namespace std;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//--------- Parse_Str_Str --------//
int Parse_Str_Str(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		if(buf=="")continue;
		out.push_back(buf);
		count++;
	}
	return count;
}
//--------- Parse_Str_Str (automatic) --------//
int Parse_Str_Str(string &in,vector <string> &out)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(! (www>>buf) )break;
		if(buf=="")continue;
		out.push_back(buf);
		count++;
	}
	return count;
}


//-------- load FASTA file -------//
int Read_FASTA_SEQRES(string &seqfile,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(seqfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",seqfile.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",seqfile.c_str());
			exit(-1);
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	seqres=temp;
	//return
	return (int)seqres.length();
}

//------------ load GFF format -------------//
//-> example
/*
##gff-version 3
##sequence-region Q6ZRR5 1 245
Q6ZRR5  UniProtKB       Chain   1       245     .       .       .       ID=PRO_0000285583;Note=Transmembrane protein 136
Q6ZRR5  UniProtKB       Transmembrane   1       21      .       .       .       Note=Helical;Ontology_term=ECO:0000255;evidence=ECO:0000255
Q6ZRR5  UniProtKB       Transmembrane   38      58      .       .       .       Note=Helical;Ontology_term=ECO:0000255;evidence=ECO:0000255
Q6ZRR5  UniProtKB       Transmembrane   75      95      .       .       .       Note=Helical;Ontology_term=ECO:0000255;evidence=ECO:0000255
Q6ZRR5  UniProtKB       Transmembrane   99      119     .       .       .       Note=Helical;Ontology_term=ECO:0000255;evidence=ECO:0000255
Q6ZRR5  UniProtKB       Transmembrane   162     182     .       .       .       Note=Helical;Ontology_term=ECO:0000255;evidence=ECO:0000255
Q6ZRR5  UniProtKB       Transmembrane   191     211     .       .       .       Note=Helical;Ontology_term=ECO:0000255;evidence=ECO:0000255
Q6ZRR5  UniProtKB       Domain  29      204     .       .       .       Note=TLC;Ontology_term=ECO:0000255;evidence=ECO:0000255|PROSITE-ProRule:PRU00205
Q6ZRR5  UniProtKB       Alternative sequence    1       1       .       .       .       ID=VSP_039853;Note=In isoform 3 and isoform 4. M->MTQCCFLRVHPLFFWFWSFHHRM;Ontology_term=ECO:0000303,ECO:0000303;evidence=ECO:0000303|PubMed:14702039,ECO:0000303|PubMed:15489334;Dbxref=PMID:14702039,PMID:15489334
Q6ZRR5  UniProtKB       Alternative sequence    91      210     .       .       .       ID=VSP_024871;Note=In isoform 2 and isoform 4. WCVYFQSEGALMLAHHTLSILGIIMALVLGESGTEVNAVLFGSELTNPLLQMRWFLRETGHYHSFTGDVVDFLFVALFTGVRIGVGACLLFCEMVSPTPKWFVKAGGVAMYAVSWCFMFS->C;Ontology_term=ECO:0000303,ECO:0000303;evidence=ECO:0000303|PubMed:14702039,ECO:0000303|PubMed:15489334;Dbxref=PMID:14702039,PMID:15489334
Q6ZRR5  UniProtKB       Sequence conflict       102     102     .       .       .       Note=M->I;Ontology_term=ECO:0000305;evidence=ECO:0000305
Q6ZRR5  UniProtKB       Sequence conflict       130     130     .       .       .       Note=L->P;Ontology_term=ECO:0000305;evidence=ECO:0000305
....
*/

int Load_GFF_TransMemb(string &fn, vector <pair<int,int> > &out)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fn.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",fn.c_str());
		exit(-1);
	}
	//load
	int col_size=9;
	int count=0;
	out.clear();
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf[0]=='#')continue;
		vector <string> tmp_str;
		int retv=Parse_Str_Str(buf,tmp_str,'\t');
		if(retv!=col_size)
		{
			fprintf(stderr,"retv %d not equal to col_size %d \n",
				retv,col_size);
			exit(-1);
		}
		if(tmp_str[2]!="Transmembrane")continue;
		int start=atoi(tmp_str[3].c_str())-1;
		int end=atoi(tmp_str[4].c_str())-1;
		out.push_back(pair<int,int>(start,end));
		count++;
	}
	//return
	return count;
}

//------ main process --------//
void Main_Process(string &fasta, string &gff)
{
	//-> load fasta file to get the sequence length
	string sequence;
	int seq_len=Read_FASTA_SEQRES(fasta,sequence);
	//-> load GFF file to get the range of transmembrane region
	vector <pair<int,int> > tm_range;
	int tm_num=Load_GFF_TransMemb(gff,tm_range);
	//-> assing TM label
	vector <int> label(seq_len,0);
	for(int i=0;i<tm_num;i++)
	{
		int start=tm_range[i].first;
		int end=tm_range[i].second;
		if(start>end)
		{
			fprintf(stderr,"start %d larger than end %d\n",start,end);
			exit(-1);
		}
		if(start<0 || end>=seq_len)
		{
			fprintf(stderr,"start %d or end %d over-range [0, %d] \n",
				start,end,seq_len);
			exit(-1);
		}
		//--- assign ---//
		for(int k=start;k<=end;k++)label[k]=1;
	}
	//-> output label
	string name;
	getBaseName(fasta,name,'/','.');
	printf(">%s\n",name.c_str());
	printf("%s\n",sequence.c_str());
	for(int i=0;i<seq_len;i++)printf("%d",label[i]);
	printf("\n");
}


//------------ main -------------//
int main(int argc, char** argv)
{
	//---- TM2_Trans ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"GFF_TransMemb <seq_file> <gff_file> \n");
			exit(-1);
		}
		string seq_file=argv[1];
		string gff_file=argv[2];
		//process
		Main_Process(seq_file,gff_file);
		//exit
		exit(0);
	}
}

