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
	char acronym[20];
	char master[100];
	char home[100];
	char script[100];
	char data[100];
	char src[50];
	FILE *f, *g, *h;
	int p, grp, R, batch, batch_size, ds, where;

//	p = 1000;
//	grp = 30;
	p = 10000;
	grp = 200;
	R = 100;
	batch_size = 2;

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
	sprintf(data,"%s/datasets",master);

	strcpy(acronym,"DG");
	sprintf(home,"%s/DG",master);
	sprintf(script,"%s/script",home);
	strcpy(src,"DG.R");

	sprintf(fname,"%s%d",acronym,p);
	h = fopen(fname,"w");
	chmod(fname,0755);

	for ( ds=0 ; ds<5 ; ds=ds++ )
	{
		sprintf(fname,"%s/%s%d_%d",script,acronym,p,ds+1);
		sprintf(line,"%s\n",fname);
		fputs(line,h);

		g = fopen(fname,"w");
		chmod(fname,0775);

		for ( batch=0 ; batch<R ; batch+=batch_size )
		{
			sprintf(fname,"%s/%s%d_%d_%03d",script,acronym,p,ds+1,batch+1);

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
			sprintf(line,"head = \"%s/p%d_%d\"\n",data,p,ds+1);
			fputs(line,f);
			sprintf(line,"R = %d\n",batch_size);
			fputs(line,f);
			sprintf(line,"seed = %d\n",(ds%2+1)*100);
			fputs(line,f);
			sprintf(line,"p = %d\n",p);
			fputs(line,f);
			sprintf(line,"g = %d\n",grp);
			fputs(line,f);
			fputs("gsm = 30\n",f);
			fputs("n = 80\n",f);
			fputs("sigma2 = 2\n",f);
			sprintf(line,"batch = %d\n",batch);
			fputs(line,f);
			if ( ds==0 )
				fputs("DG_batch(R,head,seed,p,g,gsm,n,n,n,sigma2,Gmode=0,batch=batch)\n",f);
			if ( ds==1 )
				fputs("DG_batch(R,head,seed,p,g,gsm,n,n,n,sigma2,ediu=0,Gmode=0,batch=batch)\n",f);
			if ( ds==2 )
				fputs("DG_batch(R,head,seed,p,g,gsm,n,n,n,sigma2,Gmode=1,batch=batch)\n",f);
			if ( ds==3 )
				fputs("DG_batch(R,head,seed,p,g,gsm,n,n,n,sigma2,ediu=0,Gmode=1,batch=batch)\n",f);
			if ( ds==4 )
				fputs("DG_batch(R,head,seed,p,g,gsm,n,n,n,sigma2,Gmode=2,batch=batch)\n",f);
			fclose(f);
		}

		fclose(g);
	}
}



