#ifndef types_h
#define types_h

//---- types -------------------------------------------------------------------

typedef char           Char_t;      //Signed Character 1 byte
typedef unsigned char  UChar_t;     //Unsigned Character 1 byte
typedef short          Short_t;     //Signed Short integer 2 bytes
typedef unsigned short UShort_t;    //Unsigned Short integer 2 bytes
#ifdef R__INT16
typedef long           Int_t;       //Signed integer 4 bytes
typedef unsigned long  UInt_t;      //Unsigned integer 4 bytes
#else
typedef int            Int_t;       //Signed integer 4 bytes
typedef unsigned int   UInt_t;      //Unsigned integer 4 bytes
#endif
#ifdef R__B64
typedef int            Seek_t;      //File pointer
typedef long           Long_t;      //Signed long integer 4 bytes
typedef unsigned long  ULong_t;     //Unsigned long integer 4 bytes
#else
typedef int            Seek_t;      //File pointer
typedef long           Long_t;      //Signed long integer 8 bytes
typedef unsigned long  ULong_t;     //Unsigned long integer 8 bytes
#endif
typedef float          Float_t;     //Float 4 bytes
typedef double         Double_t;    //Float 8 bytes
typedef char           Text_t;      //General string
// typedef unsigned char  Bool_t;      //Boolean (0=false, 1=true)
typedef unsigned char  Byte_t;      //Byte (8 bits)
typedef short          Version_t;   //Class version identifier
typedef const char     Option_t;    //Option string
typedef int            Ssiz_t;      //String size
typedef float          Real_t;      //TVector and TMatrix element type

typedef void         (*VoidFuncPtr_t)();  //pointer to void function
typedef bool           Bool_t;      //Boolean (0=false, 1=true) (bool)

//---- constants ---------------------------------------------------------------

#ifndef NULL
#define NULL 0
#endif

// const Bool_t kTRUE   = 1;
// const Bool_t kFALSE  = 0;

// const Int_t  kMaxInt      = 2147483647;
// const Int_t  kMaxShort    = 32767;
const UInt_t    kMaxUInt     = UInt_t(~0);
const Int_t     kMaxInt      = Int_t(kMaxUInt >> 1);
const Int_t     kMaxUShort   = 65534;
const Int_t     kMaxShort    = kMaxUShort >> 1;

//const size_t kBitsPerByte = 8;
//const Ssiz_t kNPOS        = ~(Ssiz_t)0;

#endif
