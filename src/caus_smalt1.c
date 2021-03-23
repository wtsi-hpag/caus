/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2017  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of CAUS Chromosome Assignment Using Synteny pipeline  *
 *                                                                          *
 *  Scaff10x is a free software: you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/


#include <math.h>
#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 50000 
#define Max_N_NameBase 60
#define Max_N_Pair 100
static char **S_Name,**R_Name,**R_Name2,**T_Name,**cellname;
static int *insert,*hit_mask,*hit_rcdex,*hit_locus1,*hit_locus2,*readlength,*superlength;
static int *hit_score,*pair_loend,*pair_insert,*pair_insert2;
static int *hit_read1,*hit_read2,*PairIndex,*PairIndex1,*hit_length,*cell2tag,*refs2tag;
static int *indel2read,*indel2copy,*indel2ctgs,*insert_size;
static float *hit_identy;

/* SSAS default parameters   */
static int IMOD=0;
static int n_type=0;
static int mapNumber=2;
static int file_flag=2;
static int tiles_flag=0;
static int edge_set=2;
static int edge_flag=0;
static int nContig=0;
static int max_len = 100000;
typedef struct
{
       int foffset;
       int fsindex;
} SIO;

fasta *expt;

static char rc_char[500000];
static char rc_sub[5000];

int ReverseComplement(int seqdex)
{
        int i,len;
        char *tp,*dp;
        fasta *seqp;

        seqp=expt+seqdex;
        len=seqp->length;
        memset(rc_sub,'\0',5000);
        dp=rc_sub;      
        tp = seqp->data+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int Reverse_Complement_Contig(char c_array[],int num_len)
{
        int i,len;
        char *tp,*dp;

        len=num_len;
        dp=rc_char;
        tp = c_array+len;
        for(i=len;--i>=0;)
        {
                int tmp = *--tp;
                if     (tmp == 't') *dp++ = 'a';
                else if(tmp == 'g') *dp++ = 'c';
                else if(tmp == 'c') *dp++ = 'g';
                else if(tmp == 'a') *dp++ = 't';
                else                *dp++ = tmp;
        }
        return(0);
}


int main(int argc, char **argv)
{
    FILE *namef,*namef2;
    int i,j,nSeq,args;
    int n_contig,n_reads,n_readsMaxctg,nseq;
    fasta *seq;
    void decodeReadpair(int nSeq);
    void HashFasta_Head(int i, int nSeq);
    void HashFasta_Table(int i, int nSeq);
    void Search_SM(fasta *seq,int nSeq);
    void Assemble_SM(int arr,int brr);
    void Readname_match(fasta *seq,char **argv,int args,int nSeq,int nRead);
    void Indel_Process(char **argv,int args,int nSeq);
    void Memory_Allocate(int arr);
    char line[2000]={0},tempc1[60],cc[60],RC[10],readname[Max_N_NameBase],*st,*ed;
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    void Read_Pairs(char **argv,int args,fasta *seq,int nSeq);

    seq=NULL;
    fflush(stdout);
    system("ps aux | grep long_indel; date");
    if(argc < 2)
    {
      printf("Usage: %s [-edge 2000] <cross_genome_ouput file>\n",argv[0]);

      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-type"))
       {
         sscanf(argv[++i],"%d",&n_type); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-edge"))
       {
         sscanf(argv[++i],"%d",&edge_set);
         edge_flag=1;
         args=args+2;
       }
       else if(!strcmp(argv[i],"-tile"))
       {
         sscanf(argv[++i],"%d",&tiles_flag);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-map"))
       {
         sscanf(argv[++i],"%d",&mapNumber);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-max"))
       {
         sscanf(argv[++i],"%d",&max_len);
         args=args+2;
       }
       else if(!strcmp(argv[i],"-file"))
       {
         sscanf(argv[++i],"%d",&file_flag);
         args=args+2;
       }
    }

    nseq=0;
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: args \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nseq++;
    }
    fclose(namef); 
   
/*
    nRead=0;
    if((namef = fopen(argv[args+1],"r")) == NULL)
    {
      printf("ERROR main:: args+1 \n");
      exit(1);
    }
    while(!feof(namef))
    {
      fgets(line,2000,namef);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef);   */ 

    if((insert = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - insert\n");
      exit(1);
    }
    if((hit_mask = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_score\n");
      exit(1);
    }
    if((hit_rcdex = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_rcdex\n");
      exit(1);
    }
    if((hit_read1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_read1\n");
      exit(1);
    }
    if((hit_read2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_read2\n");
      exit(1);
    }
    if((hit_locus1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus1\n");
      exit(1);
    }
    if((hit_locus2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_locus2\n");
      exit(1);
    }
    if((hit_length = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_length\n");
      exit(1);
    }
    if((readlength = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - readlength\n");
      exit(1);
    }
    if((superlength = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - superlength\n");
      exit(1);
    }
    if((hit_score = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - hit_score\n");
      exit(1);
    }
    if((pair_loend = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - pair_loend\n");
      exit(1);
    }
    if((pair_insert = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - pair_insert\n");
      exit(1);
    }
    if((pair_insert2 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - pair_insert2\n");
      exit(1);
    }
    if((PairIndex = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - PairIndex\n");
      exit(1);
    }
    if((PairIndex1 = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - PairIndex1\n");
      exit(1);
    }
    if((hit_identy = (float *)calloc(nseq,sizeof(float))) == NULL)
    {
      printf("fmate: calloc - hit_identy\n");
      exit(1);
    }
    if((indel2read = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - indel2read\n");
      exit(1);
    }
    if((indel2ctgs = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - indel2ctgs\n");
      exit(1);
    }
    if((indel2copy = (int *)calloc(nseq,sizeof(int))) == NULL)
    {
      printf("fmate: calloc - indel2copy\n");
      exit(1);
    }

    nSeq=nseq;
    R_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
    S_Name=cmatrix(0,nseq+10,0,Max_N_NameBase);
    T_Name=cmatrix(0,nseq+10,0,6);
    n_readsMaxctg=0;
    n_contig=0;
    n_reads=0;

    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    if((namef2 = fopen(argv[args+1],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the alignment files         */
    printf("www: %d %s\n",nseq,argv[args]);
    i=0;
    while(fscanf(namef,"%s %d %s %s %d %d %d %d %s %d %f %d %d",tempc1,&hit_score[i],readname,S_Name[i],&hit_read1[i],&hit_read2[i],&hit_locus1[i],&hit_locus2[i],RC,&hit_length[i],&hit_identy[i],&readlength[i],&superlength[i])!=EOF)
    {
        int idt;
        st = readname;
        ed= strrchr(readname,'_');
	strncpy(R_Name[i],readname,ed-st);
        memset(cc,'\0',60);
        for(j=10;j<19;j++)
           cc[j-10] = readname[j]; 
        idt = atoi(cc);
        if(RC[0] == 'F')
        {
          hit_read1[i] = hit_read1[i] + idt;
          hit_read2[i] = hit_read2[i] + idt;
          hit_rcdex[i]=0;
        }
        else
	{
          int hit1 = hit_locus1[i];
	  int hit2 = hit_locus2[i];
          hit_locus1[i] = hit2;
          hit_locus2[i] = hit1;
          hit1 = hit_read1[i];
          hit2 = hit_read2[i];
          hit_read1[i] = hit2 + idt;
          hit_read2[i] = hit1 + idt;

//	  hit_locus1[i] = superlength[i] - hit1 +1;
//	  hit_locus2[i] = superlength[i] - hit2 +1;
          hit_rcdex[i]=1;
	}
	hit_mask[i] = 0;
	strncpy(T_Name[i],S_Name[i],4);
    fprintf(namef2,"ALIGNMENT %s %d %s %d %d %d %d %s %d %f %d %d\n",R_Name[i],idt,S_Name[i],hit_read1[i],hit_read2[i],hit_locus1[i],hit_locus2[i],RC,hit_length[i],hit_identy[i],readlength[i],superlength[i]);
        i++;
    }
    fclose(namef);
    fclose(namef2);

    printf("Job finished for %d reads!\n",nSeq);
    return EXIT_SUCCESS;

}
/* end of the main */

/*   subroutine to sort out read pairs    */
/* =============================== */
void Indel_Process(char **argv,int args,int nSeq)
/* =============================== */
{
     int i,j,k,m,n;
     int num_hits,hit_ray[5000];
     int stopflag,*readIndex,*readIndex2;
     int offset;
     void ArraySort_Mix(int n, long *arr, int *brr);
     char **DBname;
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     void ArraySort_String(int n,char **Pair_Name,int *brr);
          
     if((readIndex= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - readIndex\n");
       exit(1);
     }
     if((readIndex2= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - readIndex\n");
       exit(1);
     }
     DBname=cmatrix(0,nSeq,0,Max_N_NameBase);
     for(i=0;i<nSeq;i++)
     {
        strncpy(DBname[i],S_Name[i],5);
        readIndex[i] = i;
     }
     num_hits =0;
     k = 0;
     offset = 0;
     for(i=0;i<(nSeq-1);i++)
     {
        stopflag=0;
        j=i+1;
        while((j<nSeq)&&(stopflag==0))
        {
          if(strcmp(R_Name[i],R_Name[j])==0)
          {
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=4) 
        {
          int kk,stopflag2 = 0;
	  kk = j-i;
	  for(n=i;n<j;n++)
	  {
	     m = n+1;
	     stopflag2 = 0;
             while((m<j)&&(stopflag2==0))
	     {
               if((strcmp(DBname[n],DBname[m])==0)&&(m<(j)))
               {
                 m++;
               }
               else
                 stopflag2=1;
	     }
	     if(((m-n) < kk)&&((m-n)>=2))
	     {
               int n_block = 0;
               int g,h,stopflag3 = 0;
//                   printf("name1: %s %s %d %d\n",DBname[n],R_Name[n],m-n,m);
	       for(g=i;g<j;g++)
	       {
	          h = g+1;
	          stopflag3 = 0;
                  while((h<j)&&(stopflag3==0))
	          {
                    if((strcmp(DBname[g],DBname[h])==0)&&(h<(j)))
                    {
                      h++;
                    }
                    else
                      stopflag3=1;
	          }
                  hit_ray[n_block] = h-g;
                  n_block++;
                  g = h - 1;
               }
               if(n_block==2)
               {
//                   printf("name2: %s %s %d %d\n",DBname[n],R_Name[n],m-n,m);
                 if((hit_ray[0]>=edge_set)&&(hit_ray[1]>=edge_set))
                 {
//                   printf("translocation: %s %s %d\n",R_Name[n],S_Name[n],kk);
                   if((S_Name[i][0]!= '*')&&(S_Name[j-1][0]!='*'))
                   {
                     printf("translocation: %s %s %d %d\n",S_Name[i],S_Name[j-1],hit_locus2[i+hit_ray[0]-1],hit_locus1[j-hit_ray[1]]);
                     for(g=i;g<j;g++)
                        printf("contigs: %s %s %d %d %d %d %d %d %d\n",R_Name[g],S_Name[g],hit_read1[g],hit_read2[g],hit_locus1[g],hit_locus2[g],hit_rcdex[g],readlength[g],superlength[g]);
                     printf("\n");
                   }
                 }
               }
	       n = j;
	     }
             else
	       n = m-1;
          }
        }
        else
        {
//          printf("www: %s %d %d\n",R_Name[i],hit_read2[i],superlength[i]);
        }
	num_hits = j-i;
	offset = offset+num_hits;
        i=j-1;
     }
     exit(1);
}


/*   subroutine to sort out read pairs    */
/* =============================== */
void Readname_match(fasta *seq, char **argv,int args,int nSeq,int nRead)
/* =============================== */
{
     FILE *namef;
     int i,j,k=0,n_reads,i_reads;
     int *readIndex;
     int stopflag,idd = 0;
     char **DBname,tagname[Max_N_NameBase],*st,*ed;
     void ArraySort_String(int n,char **Pair_Name,int *brr);
     char **cmatrix(long nrl,long nrh,long ncl,long nch);
     int num_rd_find=0;
     char line[500],base[100],zero[100]={0},*ptr;

     n_reads = nSeq+nRead; 
     DBname=cmatrix(0,n_reads+1,0,Max_N_NameBase);
     R_Name2=cmatrix(0,nRead+1,0,Max_N_NameBase);
     cellname=cmatrix(0,nRead+1,0,Max_N_NameBase);
     if((readIndex= (int *)calloc(n_reads,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - readIndex\n");
       exit(1);
     }
     if((cell2tag= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - icell2tag\n");
       exit(1);
     }
     if((refs2tag= (int *)calloc(nSeq,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - icell2tag\n");
       exit(1);
     }
     if((insert_size= (int *)calloc(nRead,sizeof(int))) == NULL)
     {
       printf("ERROR Memory_Allocate: calloc - insert_size\n");
       exit(1);
     }
     memset(cell2tag,-1,nSeq*4);
     memset(refs2tag,-1,nSeq*4);
     if((namef = fopen(argv[args+1],"r")) == NULL)
     {
       printf("ERROR Memory_Allocate:: reads group file \n");
       exit(1);
     }
     
     i_reads = 0;
     while(!feof(namef))
     {
       fgets(line,500,namef);
       if(feof(namef)) break;
       i = 0;
       for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
       {
          if(i==0)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            strcpy(R_Name2[i_reads],base);
          }
          else if(i==2)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            insert_size[i_reads]=atoi(ptr);
          }
          else if(i==3)
          {
            strcpy(base,zero);
            strcat(base,ptr);
            strcpy(tagname,base);
            st = tagname;
            ed = strchr(tagname,'.');
//            ed = strrchr(tagname,'_');
            if(ed==NULL)
              strcpy(cellname[i_reads],tagname);
            else
              strncpy(cellname[i_reads],tagname,ed-st);
            i_reads++;
            ptr = NULL;
          }
       }
     }

/*   find out the read name match   */
     for(j=0;j<nSeq;j++)
     {
        strcpy(DBname[j],R_Name[j]);
        readIndex[j]=j;
     }
     for(j=0;j<nRead;j++)
     {
        strcpy(DBname[j+nSeq],R_Name2[j]);
        readIndex[j+nSeq]=j+nSeq;
     }
     n_reads=nSeq+nRead;
     ArraySort_String(n_reads,DBname,readIndex);

     num_rd_find=0;
     for(i=0;i<n_reads;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        if(readIndex[i]>=nSeq)
          idd = readIndex[i];
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            if(readIndex[j]>=nSeq)
              idd = readIndex[j];
//            num_rd_find++;
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2)
        {
          for(k=i;k<j;k++)
          {
             if(readIndex[k]<nSeq)
             {
               cell2tag[readIndex[k]] = idd-nSeq;
//         printf("name: %d %s\n",k,DBname[k]);
               num_rd_find++;
             }
          }
        }
        i=j-1;
     }
     printf("reads found: %d %d %d\n",nSeq,num_rd_find,nRead);

/*   find out the read name match   */
     for(j=0;j<nSeq;j++)
     {
        strcpy(DBname[j],S_Name[j]);
        readIndex[j]=j;
     }
     nRead = nContig;
     for(j=0;j<nRead;j++)
     {
        strcpy(DBname[j+nSeq],(seq+j)->name);
        readIndex[j+nSeq]=j+nSeq;
     }
     n_reads=nSeq+nRead;
     ArraySort_String(n_reads,DBname,readIndex);

     num_rd_find=0;
     idd = 0;
     for(i=0;i<n_reads;i++)
     {
/*      search reads with an index < i     */
/*      search reads with an index > i     */
        stopflag=0;
        j=i+1;
        if(readIndex[i]>=nSeq)
          idd = readIndex[i];
        while((j<n_reads)&&(stopflag==0))
        {
          if(strcmp(DBname[j],DBname[i])==0)
          {
            if(readIndex[j]>=nSeq)
              idd = readIndex[j];
//            num_rd_find++;
            j++;
          }
          else
            stopflag=1;
        }
        if((j-i)>=2)
        {
          for(k=i;k<j;k++)
          {
             if(readIndex[k]<nSeq)
             {
               refs2tag[readIndex[k]] = idd-nSeq;
//         printf("name: %d %s\n",k,DBname[k]);
               num_rd_find++;
             }
          }
        }
        i=j-1;
     }
     printf("contigs found: %d %d %d\n",nSeq,num_rd_find,nRead);
}


#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char **Pair_Name, int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char **Pair_Name, int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

