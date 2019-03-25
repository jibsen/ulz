/*

ULZ.HPP - An ultra-fast LZ77 data compression library

Written and placed in the public domain by Ilya Muravyov

*/

#ifndef ULZ_HPP_INCLUDED
#define ULZ_HPP_INCLUDED

class CULZ
{
public:
  typedef unsigned char U8;
  typedef unsigned short U16;
  typedef unsigned int U32;

  static const int EXCESS=16;

  static const int WINDOW_BITS=17;
  static const int WINDOW_SIZE=1<<WINDOW_BITS;
  static const int WINDOW_MASK=WINDOW_SIZE-1;

  static const int MIN_MATCH=4;

  static const int HASH_BITS=18;
  static const int HASH_SIZE=1<<HASH_BITS;
  static const int NIL=-1;

  // Hash chain

  int Head[HASH_SIZE];
  int Tail[WINDOW_SIZE];

  // Utils

  template<typename T>
  inline T Min(T a, T b)
  {
    return (a<b)?a:b;
  }

  template<typename T>
  inline T Max(T a, T b)
  {
    return (a>b)?a:b;
  }

  inline U16 UnalignedLoad16(void* p)
  {
    return *reinterpret_cast<const U16*>(p);
  }

  inline U32 UnalignedLoad32(void* p)
  {
    return *reinterpret_cast<const U32*>(p);
  }

  inline void UnalignedStore16(void* p, U16 x)
  {
    *reinterpret_cast<U16*>(p)=x;
  }

  inline void UnalignedCopy32(void* d, void* s)
  {
    *reinterpret_cast<U32*>(d)=UnalignedLoad32(s);
  }

  inline void WildCopy(U8* d, U8* s, int n)
  {
    UnalignedCopy32(d, s);
    UnalignedCopy32(d+4, s+4);

    for (int i=8; i<n; i+=8)
    {
      UnalignedCopy32(d+i, s+i);
      UnalignedCopy32(d+i+4, s+i+4);
    }
  }

  inline U32 Hash32(void* p)
  {
    return (UnalignedLoad32(p)*0x9E3779B9)>>(32-HASH_BITS);
  }

  inline void EncodeMod(U8*& p, U32 x)
  {
    while (x>=128)
    {
      x-=128;
      *p++=128+(x&127);
      x>>=7;
    }
    *p++=x;
  }

  inline U32 DecodeMod(U8*& p)
  {
    U32 x=0;
    for (int i=0; i<=28; i+=7)
    {
      const U32 c=*p++;
      x+=c<<i;
      if (c<128)
        break;
    }
    return x;
  }

  // LZ77

  int Compress(U8* in, int in_len, U8* out, int level=4)
  {
    const int max_chain=(level<8)?1<<level:WINDOW_SIZE;
    U8* op=out;

    for (int i=0; i<HASH_SIZE; ++i)
      Head[i]=NIL;

    int run=0;

    int p=0;
    while (p<in_len)
    {
      int best_len=MIN_MATCH-1;
      int dist=0;

      const int max_match=in_len-p;
      if (max_match>=MIN_MATCH)
      {
        int limit=Max(p-WINDOW_SIZE, NIL);
        int chain_len=max_chain;

        int s=Head[Hash32(&in[p])];
        while (s>limit)
        {
          if (in[s+best_len]==in[p+best_len]
              && UnalignedLoad32(&in[s])==UnalignedLoad32(&in[p]))
          {
            int len=MIN_MATCH;
            while (len<max_match && in[s+len]==in[p+len])
              ++len;

            if (len>best_len)
            {
              best_len=len;
              dist=p-s;

              if (len==max_match)
                break;
            }
          }

          if (!--chain_len)
            break;

          s=Tail[s&WINDOW_MASK];
        }

        if (best_len==MIN_MATCH && run>=(7+128))
          best_len=0;

        if (level==9 && best_len>=MIN_MATCH && best_len<max_match)
        {

          for (int i=1; i<=2 && best_len; ++i) // 2-byte lookahead
          {
            const int target_len=best_len+i;
            const int j=p+i;

            limit=Max(j-WINDOW_SIZE, NIL);
            chain_len=max_chain;

            s=Head[Hash32(&in[j])];
            while (s>limit)
            {
              if (in[s+best_len]==in[j+best_len]
                  && UnalignedLoad32(&in[s])==UnalignedLoad32(&in[j]))
              {
                int len=MIN_MATCH;
                while (len<target_len && in[s+len]==in[j+len])
                  ++len;

                if (len==target_len)
                {
                  best_len=0;
                  break;
                }
              }

              if (!--chain_len)
                break;

              s=Tail[s&WINDOW_MASK];
            }
          }
        }

      }

      if (best_len>=MIN_MATCH)
      {
        const int len=best_len-MIN_MATCH;
        const int tmp=((dist>>12)&16)+Min(len, 15);

        if (run)
        {
          if (run>=7)
          {
            *op++=(7<<5)+tmp;
            EncodeMod(op, run-7);
          }
          else
            *op++=(run<<5)+tmp;

          WildCopy(op, &in[p-run], run);
          op+=run;

          run=0;
        }
        else
          *op++=tmp;

        if (len>=15)
          EncodeMod(op, len-15);

        UnalignedStore16(op, dist);
        op+=2;

        while (best_len--)
        {
          const U32 h=Hash32(&in[p]);
          Tail[p&WINDOW_MASK]=Head[h];
          Head[h]=p++;
        }
      }
      else
      {
        ++run;

        const U32 h=Hash32(&in[p]);
        Tail[p&WINDOW_MASK]=Head[h];
        Head[h]=p++;
      }
    }

    if (run)
    {
      if (run>=7)
      {
        *op++=7<<5;
        EncodeMod(op, run-7);
      }
      else
        *op++=run<<5;

      WildCopy(op, &in[p-run], run);
      op+=run;

      //run=0;
    }

    return op-out;
  }

  int Decompress(U8* in, int in_len, U8* out, int out_len)
  {
    U8* op=out;
    U8* ip=in;
    const U8* ip_end=ip+in_len;
    const U8* op_end=op+out_len;

    while (ip<ip_end)
    {
      const U32 tag=*ip++;
      if (tag>=32)
      {
        U32 run=tag>>5;
        if (run==7)
          run+=DecodeMod(ip);
        if ((op_end-op)<run || (ip_end-ip)<run) // Overrun check
          return -1;

        WildCopy(op, ip, run);
        op+=run;
        ip+=run;
        if (ip>=ip_end)
          break;
      }

      U32 len=(tag&15)+MIN_MATCH;
      if (len==(15+MIN_MATCH))
        len+=DecodeMod(ip);
      if ((op_end-op)<len) // Overrun check
        return -1;

      const U32 dist=((tag&16)<<12)+UnalignedLoad16(ip);
      ip+=2;
      if ((op-out)<dist) // Range check
        return -1;
      U8* cp=op-dist;

      if (dist>=4)
      {
        WildCopy(op, cp, len);
        op+=len;
      }
      else
      {
        while (len--)
          *op++=*cp++;
      }
    }

    return (ip==ip_end)?op-out:-1;
  }
};

#endif // ULZ_HPP_INCLUDED
