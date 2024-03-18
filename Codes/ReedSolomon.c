#include <stdio.h>
#include <time.h>
// 14 16 18 20
#define mm 4
/* GF(2**3) */
#define nn 14
/* nn=2**mm -1 length of codeword */
#define tt 1
/* number of errors that can be corrected */
#define kk 12
/* kk = nn-2*tt */
int nr_poly[mm + 1] = {1, 1, 0, 0,1};
/* specify irreducible polynomial coeffts */
int al[nn + 1], index_of[nn + 1], gg[nn - kk + 1];
int recd[nn], data[kk], bb[nn - kk];
char MAC[kk];

float total_time = 0;

int convert_to_num(char c){
    int num = c-'0';
    if(num<10){
        return num;
    }
    else if(num==17 || num==49)
    return 10;
    else if(num==18 || num==50)
    return 11;
    else if(num==19 || num==51)
    return 12;
    else if(num==20 || num==52)
    return 13;
    else if(num==21 || num==53)
    return 14;
    else if(num==22 || num==54)
    return 15;
    
}

char convert_to_char(int n){
    if(n==0)return '0';
    if(n==1)return '1';
    else if(n==2)return '2';
    else if(n==3)return '3';
    else if(n==4)return '4';
    else if(n==5)return '5';
    else if(n==6)return '6';
    else if(n==7)return '7';
    else if(n==8)return '8';
    else if(n==9)return '9';
    else if(n==10)return 'A';
    else if(n==11){
        return 'B';
    }
    else if(n==12)return 'C';
    else if(n==13)return 'D';
    else if(n==14)return 'E';
    else if(n==15)return 'F';
}


void generate_gf() {
    register int i, mask;
    printf("void generate_gf()");
    mask = 1;
    al[mm] = 0;
    for (i = 0; i < mm; i++) {
        al[i] = mask;
        index_of[al[i]] = i;
        if (nr_poly[i] != 0)
            al[mm] ^= mask;
        mask <<= 1;
    }
    index_of[al[mm]] = mm;
    mask >>= 1;
    for (i = mm + 1; i < nn; i++) {
        if (al[i - 1] >= mask)
            al[i] = al[mm] ^ ((al[i - 1] ^ mask) << 1);
        else al[i] = al[i - 1] << 1;
        index_of[al[i]] = i;
    }
    index_of[0] = -1;
}

void gen_poly() {
    register int i, j;
    printf("\n void gen_poly()");
    gg[0] = 2;
    /* primitive element alpha = 2 for GF(2**mm) */
    gg[1] = 1;
    /* g(x) = (X+alpha) initially */
    for (i = 2; i <= nn - kk; i++) {
        gg[i] = 1;
        for (j = i - 1; j > 0; j--)
            if (gg[j] != 0) gg[j] = gg[j - 1] ^ al[(index_of[gg[j]] + i) % nn];
            else gg[j] = gg[j - 1];
        gg[0] = al[(index_of[gg[0]] + i) % nn];
        /* gg[0] can never be zero */
    }
    /* convert gg[] to index form for quicker encoding */
    for (i = 0; i <= nn - kk; i++) {
        printf("\n gp gg[%2d ]=%3d ", i, gg[i]); /*display generator polynomial array */
        gg[i] = index_of[gg[i]];
    }
}

void encode_rs() {
    register int i, j;
    printf("\n void encode_rs()");
    int feedback;
    for (i = 0; i < nn - kk; i++) bb[i] = 0;
    for (i = kk - 1; i >= 0; i--) {
        feedback = index_of[data[i] ^ bb[nn - kk - 1]];
        if (feedback != -1) {
            for (j = nn - kk - 1; j > 0; j--)
                if (gg[j] != -1)
                    bb[j] = bb[j - 1] ^ al[(gg[j] + feedback) % nn];
                else
                    bb[j] = bb[j - 1];
            bb[0] = al[(gg[0] + feedback) % nn];
        } else {
            for (j = nn - kk - 1; j > 0; j--)
                bb[j] = bb[j - 1];
            bb[0] = 0;
        }
    }
}
void decode_rs()
 {
   register int i,j,u,q ;
   int elp[nn-kk+2][nn-kk], d[nn-kk+2], l[nn-kk+2], u_lu[nn-kk+2], s[nn-kk+1] ;
   int count=0, syn_error=0, root[tt], loc[tt], z[tt+1], err[nn], reg[tt+1] ;

/* first form the syndromes */
   for (i=1; i<=nn-kk; i++)
    { s[i] = 0 ;
      for (j=0; j<nn; j++)
        if (recd[j]!=-1)
          s[i] ^= al[(recd[j]+i*j)%nn] ;      /* recd[j] in index form */
/* convert syndrome from polynomial form to index form  */
      if (s[i]!=0)  syn_error=1 ;        /* set flag if non-zero syndrome => error */
      s[i] = index_of[s[i]] ;
    } ;

   if (syn_error)       /* if errors, try and correct */
    {
/* initialise table entries */
      d[0] = 0 ;           /* index form */
      d[1] = s[1] ;        /* index form */
      elp[0][0] = 0 ;      /* index form */
      elp[1][0] = 1 ;      /* polynomial form */
      for (i=1; i<nn-kk; i++)
        { elp[0][i] = -1 ;   /* index form */
          elp[1][i] = 0 ;   /* polynomial form */
        }
      l[0] = 0 ;
      l[1] = 0 ;
      u_lu[0] = -1 ;
      u_lu[1] = 0 ;
      u = 0 ;

      do
      {
        u++ ;
        if (d[u]==-1)
          { l[u+1] = l[u] ;
            for (i=0; i<=l[u]; i++)
             {  elp[u+1][i] = elp[u][i] ;
                elp[u][i] = index_of[elp[u][i]] ;
             }
          }
        else
/* search for words with greatest u_lu[q] for which d[q]!=0 */
          { q = u-1 ;
            while ((d[q]==-1) && (q>0)) q-- ;
/* have found first non-zero d[q]  */
            if (q>0)
             { j=q ;
               do
               { j-- ;
                 if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                   q = j ;
               }while (j>0) ;
             } ;

/* have now found q such that d[u]!=0 and u_lu[q] is maximum */
/* store degree of new elp polynomial */
            if (l[u]>l[q]+u-q)  l[u+1] = l[u] ;
            else  l[u+1] = l[q]+u-q ;

/* form new elp(x) */
            for (i=0; i<nn-kk; i++)    elp[u+1][i] = 0 ;
            for (i=0; i<=l[q]; i++)
              if (elp[q][i]!=-1)
                elp[u+1][i+u-q] = al[(d[u]+nn-d[q]+elp[q][i])%nn] ;
            for (i=0; i<=l[u]; i++)
              { elp[u+1][i] ^= elp[u][i] ;
                elp[u][i] = index_of[elp[u][i]] ;
              }
          }
        u_lu[u+1] = u-l[u+1] ;

/* form (u+1)th discrepancy */
        if (u<nn-kk)    /* no discrepancy computed on last iteration */
          {
            if (s[u+1]!=-1)
                   d[u+1] = al[s[u+1]] ;
            else
              d[u+1] = 0 ;
            for (i=1; i<=l[u+1]; i++)
              if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                d[u+1] ^= al[(s[u+1-i]+index_of[elp[u+1][i]])%nn] ;
            d[u+1] = index_of[d[u+1]] ;    /* put d[u+1] into index form */
          }
      } while ((u<nn-kk) && (l[u+1]<=tt)) ;

      u++ ;
      if (l[u]<=tt)         /* can correct error */
       {
/* put elp into index form */
         for (i=0; i<=l[u]; i++)   elp[u][i] = index_of[elp[u][i]] ;

/* find roots of the error location polynomial */
         for (i=1; i<=l[u]; i++)
           reg[i] = elp[u][i] ;
         count = 0 ;
         for (i=1; i<=nn; i++)
          {  q = 1 ;
             for (j=1; j<=l[u]; j++)
              if (reg[j]!=-1)
                { reg[j] = (reg[j]+j)%nn ;
                  q ^= al[reg[j]] ;
                } ;
             if (!q)        /* store root and error location number indices */
              { root[count] = i;
                loc[count] = nn-i ;
                count++ ;
              };
          } ;
         if (count==l[u])    /* no. roots = degree of elp hence <= tt errors */
          {
/* form polynomial z(x) */
           for (i=1; i<=l[u]; i++)        /* Z[0] = 1 always - do not need */
            { if ((s[i]!=-1) && (elp[u][i]!=-1))
                 z[i] = al[s[i]] ^ al[elp[u][i]] ;
              else if ((s[i]!=-1) && (elp[u][i]==-1))
                      z[i] = al[s[i]] ;
                   else if ((s[i]==-1) && (elp[u][i]!=-1))
                          z[i] = al[elp[u][i]] ;
                        else
                          z[i] = 0 ;
              for (j=1; j<i; j++)
                if ((s[j]!=-1) && (elp[u][i-j]!=-1))
                   z[i] ^= al[(elp[u][i-j] + s[j])%nn] ;
              z[i] = index_of[z[i]] ;         /* put into index form */
            } ;

  /* evaluate errors at locations given by error location numbers loc[i] */
           for (i=0; i<nn; i++)
             { err[i] = 0 ;
               if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                 recd[i] = al[recd[i]] ;
               else  recd[i] = 0 ;
             }
           for (i=0; i<l[u]; i++)    /* compute numerator of error term first */
            { err[loc[i]] = 1;       /* accounts for z[0] */
              for (j=1; j<=l[u]; j++)
                if (z[j]!=-1)
                  err[loc[i]] ^= al[(z[j]+j*root[i])%nn] ;
              if (err[loc[i]]!=0)
               { err[loc[i]] = index_of[err[loc[i]]] ;
                 q = 0 ;     /* form denominator of error term */
                 for (j=0; j<l[u]; j++)
                   if (j!=i)
                     q += index_of[1^al[(loc[j]+root[i])%nn]] ;
                 q = q % nn ;
                 err[loc[i]] = al[(err[loc[i]]-q+nn)%nn] ;
                 recd[loc[i]] ^= err[loc[i]] ;  /*recd[i] must be in polynomial form */
               }
            }
          }
         else    /* no. roots != degree of elp => >tt errors and cannot solve */
           for (i=0; i<nn; i++)        /* could return error flag if desired */
               if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                 recd[i] = al[recd[i]] ;
               else  recd[i] = 0 ;     /* just output received codeword as is */
       }
     else         /* elp has degree has degree >tt hence cannot solve */
       for (i=0; i<nn; i++)       /* could return error flag if desired */
          if (recd[i]!=-1)        /* convert recd[] to polynomial form */
            recd[i] = al[recd[i]] ;
          else  recd[i] = 0 ;     /* just output received codeword as is */
    }
   else       /* no non-zero syndromes => no errors: output received codeword */
    for (i=0; i<nn; i++)
       if (recd[i]!=-1)        /* convert recd[] to polynomial form */
         recd[i] = al[recd[i]] ;
       else  recd[i] = 0 ;
 }


int main() {
    register int i, j;
    clock_t start, end;
    generate_gf();
    gen_poly();
    double total_time = 0;
    int testCases;
    printf("Enter number of MAC addresses to encode and decode : ");
    scanf("%d", &testCases);

    for (j = 0; j < testCases; j++) {
        printf("Enter MAC address: ");
        for (i = 0; i < kk; i++) {
            scanf(" %c", &MAC[i]);
        }
        for (i = 0; i < kk; i++) {
            data[i] = convert_to_num(MAC[i]);
        }
        encode_rs();
        for (i = 0; i < nn - kk; i++) recd[i] = bb[i];
        for (i = 0; i < kk; i++) recd[i + nn - kk] = data[i];

        start = clock();
        decode_rs();
        end = clock();

        double time_taken = ((double)(end - start)) / (CLOCKS_PER_SEC / 1000000);
        total_time += time_taken;
    }

    printf("Total time taken to decode %d MAC addresses: %f microseconds\n", testCases, total_time);

    return 0;
}
