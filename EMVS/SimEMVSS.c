#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define HPC 0
#define LPC 1
#define Emory 2

int main()
{
	char fname[100];
	char Rname[100];
	char line[1024];
	char method[20];
	char acronym[20];
	char master[100];
	char home[100];
	char script[100];
	char data[100];
	char src[50];
	char vname[50];
	FILE *f, *g, *h, *m;
	int p, R, s, batch_size, batch, where;

	p = 1000;
	R = 500;
	batch_size = 1;

	strcpy(method,"EMVSS");

	if ( access("/home/cchan40",X_OK) == 0 )
	{
		where = Emory;
		strcpy(master,"/home/cchan40/project/EMSHS");
	}
	else if ( access("/project/qlonglab/changgee",X_OK) == 0 )
	{
		where = HPC;
		strcpy(master,"/home/changgee/project/EMSHS");
	}
	else
	{
		where = LPC;
		strcpy(master,"/home/changgee/project/EMSHS");
	}
	sprintf(home,"%s/EMVS",master);
	sprintf(script,"%s/Sim%d",home,p);
	strcpy(src,"SimEMVSS.R");

	sprintf(data,"%s/datasets",master);

	sprintf(fname,"%s%d",method,p);
	h = fopen(fname,"w");
	chmod(fname,0755);

	sprintf(fname,"%s%dMERGE",method,p);
	m = fopen(fname,"w");
	chmod(fname,0755);

	for ( s=0 ; s<5 ; s++ )
	{
		sprintf(acronym,"%s_%d_%d",method,p,s+1);
		sprintf(vname,"res%s",acronym);

		sprintf(fname,"%s/%s",script,acronym);
		sprintf(line,"%s\n",fname);
		fputs(line,h);

		g = fopen(fname,"w");
		chmod(fname,0755);

		for ( batch=0 ; batch<R ; batch+=batch_size )
		{
			sprintf(fname,"%s/%s_%03d",script,acronym,batch+1);
			if ( where == Emory )
				sprintf(line,"qsub -q fruit.q %s\n",fname);
			else if ( where == HPC )
				sprintf(line,"bsub -q qlonglab -e %s.e -o %s.o < %s\n",fname,fname,fname);
			else
				sprintf(line,"bsub -q cceb_normal -e %s.e -o %s.o < %s\n",fname,fname,fname);
			fputs(line,g);

			f = fopen(fname,"w");
			if ( where == LPC )
				fputs("module load R\n",f);
			else if ( where == HPC )
			{
				fputs("source /etc/profile.d/modules.sh\n",f);
				fputs("module load R-3.3.1\n",f);
			}
			sprintf(Rname,"%s.R",fname);
			sprintf(line,"R --vanilla < %s\n",Rname);
			fputs(line,f);
			fclose(f);

			f = fopen(Rname,"w");
			fputs("library(compiler)\n",f);
			fputs("enableJIT(3)\n",f);
			sprintf(line,"source(\"%s/%s\")\n",home,src);
			fputs(line,f);

			sprintf(line,"r = %d\n",batch_size);
			fputs(line,f);
			sprintf(line,"datapath = \"%s/p%d_%d\"\n",data,p,s+1);
			fputs(line,f);

			if ( p == 1000 )
				fputs("v0 = exp(seq(log(0.0004),log(0.1),length.out=20))\n",f);
			else
				fputs("v0 = exp(seq(log(0.00008),log(0.002),length.out=20))\n",f);
			fputs("v1 = 1000\n",f);
			fputs("eta = 0:1\n",f);

			sprintf(line,"if ( !file.exists(\"%s/%s_%03d\") )\n",script,vname,batch+1);
			fputs(line,f);
			fputs("{\n",f);
			sprintf(line,"  %s = SimEMVSS(r,v0,v1,eta,datapath,batch=%d)\n",vname,batch);
			fputs(line,f);

			sprintf(line,"  save(%s,file=\"%s/%s_%03d\")\n",vname,script,vname,batch+1);
			fputs(line,f);
			fputs("}\n",f);
			fclose(f);
		}
		fclose(g);

		sprintf(fname,"%s/%s_MERGE",script,acronym);
		if ( where == Emory )
			sprintf(line,"qsub -q fruit.q %s\n",fname);
		else if ( where == HPC )
			sprintf(line,"bsub -q qlonglab -e %s.e -o %s.o < %s\n",fname,fname,fname);
		else
			sprintf(line,"bsub -q cceb_normal -e %s.e -o %s.o < %s\n",fname,fname,fname);
		fputs(line,m);

		f = fopen(fname,"w");
		if ( where == LPC )
			fputs("module load R\n",f);
		else if ( where == HPC )
		{
			fputs("source /etc/profile.d/modules.sh\n",f);
			fputs("module load R-3.3.1\n",f);
		}
		sprintf(Rname,"%s.R",fname);
		sprintf(line,"R --vanilla < %s\n",Rname);
		fputs(line,f);
		fclose(f);

		g = fopen(Rname,"w");

		fputs("library(abind)\n",g);
		sprintf(line,"for ( i in 1:%d )\n",R/batch_size);
		fputs(line,g);
		fputs("{\n",g);
		sprintf(line,"  fname = sprintf(\"%s/%s_%%03d\",(i-1)*%d+1)\n",script,vname,batch_size);
		fputs(line,g);
		fputs("  load(fname)\n",g);
		fputs("  if ( i==1 )\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp = %s\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("  else\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp$r = tmp$r + %s$r\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$FNrate = abind(tmp$FNrate,%s$FNrate)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$FPrate = abind(tmp$FPrate,%s$FPrate)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$MSTE = abind(tmp$MSTE,%s$MSTE)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$MSPE = abind(tmp$MSPE,%s$MSPE)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$time = abind(tmp$time,%s$time)\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("}\n",g);
		sprintf(line,"%s = tmp\n",vname);
		fputs(line,g);
		sprintf(line,"save(%s,file=\"%s/%s\")\n",vname,home,vname);
		fputs(line,g);

		fclose(g);
	}
	fclose(h);
	fclose(m);
}



