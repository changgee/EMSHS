#include <stdio.h>
#include <string.h>

int main()
{
	char fname[100];
	char line[1024];
	char master[100];
	char method[20];
	char acronym[20];
	char home[100];
	char script[100];
	char datadir[100];
	char data[20];
	char src[50];
	char vname[50];
	FILE *f, *g;
	int i, K, s;

	K = 5;
	strcpy(data,"TCGA");
	strcpy(method,"SHW");

	strcpy(master,"/home/cchan40/project/EMSHS");
	sprintf(home,"%s/EMSH",master);
	sprintf(datadir,"%s/datasets/%s",master,data);
	sprintf(script,"%s/%s",home,data);
	strcpy(src,"DataEMSHS.R");

	for ( s=1 ; s<=6 ; s++ )
	{
		sprintf(acronym,"%sCV%d%s%02d",method,K,data,s);
		sprintf(vname,"res%s",acronym);

		sprintf(fname,"%s",acronym);
		g = fopen(fname,"w");
		chmod(fname,0755);

		for ( i=0 ; i<K ; i++ )
		{
			sprintf(line,"bsub < %s/%s%02d\n",script,acronym,i+1);
			fputs(line,g);

			sprintf(fname,"%s/%s%02d.R",script,acronym,i+1);
			f = fopen(fname,"w");
			fputs("library(compiler)\n",f);
			fputs("enableJIT(3)\n",f);
			sprintf(line,"source(\"%s/%s\")\n",home,src);
			fputs(line,f);
			sprintf(line,"load(\"%s/data\")\n",datadir);
			fputs(line,f);
			sprintf(line,"load(\"%s/fold%d\")\n",datadir,K);
			fputs(line,f);

			fputs("mu = seq(6.3,5.8,length.out=20)\n",f);
			if ( s == 1 )
			{
				fputs("nu = 0.0001\n",f);
			}
			else if ( s == 2 )
			{
				fputs("nu = 0.0003\n",f);
			}
			else if ( s == 3 )
			{
				fputs("nu = 0.001\n",f);
			}
			else if ( s == 4 )
			{
				fputs("nu = 0.003\n",f);
			}
			else if ( s == 5 )
			{
				fputs("nu = 0.01\n",f);
			}
			else if ( s == 6 )
			{
				fputs("nu = 0.03\n",f);
			}

			sprintf(line,"%s = DataEMSH_CV(y,X,mu,nu,w=w,fold=f,k=%d)\n",vname,i+1);
			fputs(line,f);

			sprintf(line,"save(%s,file=\"%s/%s%02d\")\n",vname,script,vname,i+1);
			fputs(line,f);
			fclose(f);

			sprintf(fname,"%s/%s%02d",script,acronym,i+1);
			f = fopen(fname,"w");
			fputs("module load R\n",f);
			sprintf(line,"R --vanilla < %s/%s%02d.R\n",script,acronym,i+1);
			fputs(line,f);
			fclose(f);
		}
		fclose(g);

		sprintf(fname,"%sMERGE.R",acronym);
		g = fopen(fname,"w");

		fputs("library(abind)\n",g);

		sprintf(line,"for ( i in 1:%d )\n",K);
		fputs(line,g);
		fputs("{\n",g);
		sprintf(line,"  fname = sprintf(\"%s/%s%%02d\",i)\n",script,vname);
		fputs(line,g);
		fputs("  load(fname)\n",g);
		fputs("  if ( i==1 )\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp = %s\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("  else\n",g);
		fputs("  {\n",g);
		sprintf(line,"    tmp$k = c(tmp$k,%s$k)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$beta = abind(tmp$beta,%s$beta,along=2)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$lambda = abind(tmp$lambda,%s$lambda,along=2)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$L = abind(tmp$L,%s$L,along=1)\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$SSPECV = tmp$SSPECV + %s$SSPECV\n",vname);
		fputs(line,g);
		sprintf(line,"    tmp$time = tmp$time + %s$time\n",vname);
		fputs(line,g);
		fputs("  }\n",g);
		fputs("}\n",g);
		fputs("tmp$MSPECV = tmp$SSPECV / tmp$n\n",g);
		sprintf(line,"%s = tmp\n",vname);
		fputs(line,g);
		sprintf(line,"save(%s,file=\"%s/%s\")\n",vname,home,vname);
		fputs(line,g);

		fclose(g);

		sprintf(fname,"%sMERGE",acronym);
		f = fopen(fname,"w");
		fputs("module load R\n",f);
		sprintf(line,"R --vanilla < %s/%sMERGE.R\n",home,acronym);
		fputs(line,f);
		fclose(f);
	}
}




