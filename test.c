 /* Return 1 if string is symmetrical, 0 otherwise. */
int symmetry(const char* seq) {
   register char s;
   register char e;
   const char *seq_end=seq;
   int i = 0;
   int seq_len=strlen(seq);
   int mp = seq_len/2;
   if(seq_len%2==1) {
      return 0;
   }
   seq_end+=seq_len;
   seq_end--;
   while(i<mp) {
      i++;
      s=*seq;
      e=*seq_end;
      if ((s=='A' && e!='T')
          || (s=='T' && e!='A')
          || (e=='A' && s!='T')
          || (e=='T' && s!='A')) {
         return 0;
      }
      if ((s=='C' && e!='G')
          || (s=='G' && e!='C')
          || (e=='C' && s!='G')
          || (e=='G' && s!='C')) {
         return 0;
      }
      seq++;
      seq_end--;
   }
   return 1;
}

char seq[]="DASDSA";

int res;
res=symmetry(seq);