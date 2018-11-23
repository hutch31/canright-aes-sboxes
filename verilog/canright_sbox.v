/* find either Sbox or its inverse in GF(2^8), by Canright Algorithm */
module canright_sbox ( A, encrypt, Q );
  input [7:0] A;
  input encrypt; /* 1 for Sbox, 0 for inverse Sbox */
  output [7:0] Q;

/* square in GF(2^2), using normal basis [Omega^2,Omega]
 * inverse is the same as square in GF(2^2), using any normal basis
 */
function [1:0] GF_SQ_2;
  input [1:0] A;
  begin
   GF_SQ_2 = { A[0], A[1] };
  end
endfunction

/* scale by w = Omega in GF(2^2), using normal basis [Omega^2,Omega] */
function [1:0] GF_SCLW_2;
  input [1:0] A;
  begin
    GF_SCLW_2 = { (A[1] ^ A[0]), A[1] };
  end
endfunction

/* scale by w^2 = Omega^2 in GF(2^2), using normal basis [Omega^2,Omega] */
function [1:0] GF_SCLW2_2;
  input [1:0] A;
  begin
  GF_SCLW2_2 = { A[0], (A[1] ^ A[0]) };
  end
endfunction

/* multiply in GF(2^2), shared factors, module GF_MULS_2 ( A, ab, B, cd, Q );
using normal basis [Omega^2,Omega] */
function [1:0] GF_MULS_2;
  input [1:0] A;
  input ab;
  input [1:0] B;
  input cd;
  reg abcd, p, q;
  begin
  abcd = ~(ab & cd);  /* note:~& syntax for NAND won’t compile */
  p = (~(A[1] & B[1])) ^ abcd;
  q = (~(A[0] & B[0])) ^ abcd;
  GF_MULS_2 = { p, q }; 
  end
endfunction

/* multiply & scale by N in GF(2^2), shared factors, basis [Omega^2,Omega] */
function [1:0] GF_MULS_SCL_2;
  input [1:0] A;
  input ab;
  input [1:0] B;
  input cd;
  reg t, p, q;
  begin
  t = ~(A[0] & B[0]); /* note: ~& syntax for NAND won’t compile */
  p = (~(ab & cd)) ^ t;
  q = (~(A[1] & B[1])) ^ t;
  GF_MULS_SCL_2 = { p, q };
  end
endfunction

function [3:0] GF_INV_4;
  input [3:0] A;
  reg [1:0] a, b, c,d,p,q;
  reg sa, sb, sd; /* for shared factors in multipliers */
  reg [1:0] ab, ab2, ab2N;
  begin
  
    a = A[3:2]; 
    b = A[1:0];
    sa = a[1] ^ a[0];
    sb = b[1] ^ b[0];
  /* optimize this section as shown below
  */
    ab = GF_MULS_2(a, sa, b, sb);
    ab2 = GF_SQ_2( (a ^ b));
    ab2N = GF_SCLW2_2( ab2);
    d = GF_SQ_2( (ab ^ ab2N));

//  assign c = { /* note: ~| syntax for NOR won’t compile */
//    ~(a[1] | b[1]) ^ (~(sa &  sb)) ,
//    ~(sa | sb) ^ (~(a[0] & b[0])) };


  /* end of optimization */

    sd = d[1] ^ d[0];
    p = GF_MULS_2(d, sd, b, sb);
    q = GF_MULS_2(d, sd, a, sa);
    GF_INV_4 = { p, q };
  end
endfunction

/* square & scale by nu in GF(2^4)/GF(2^2),
 * in the normal basis [alpha^8, alpha^2]:
 *
 * nu = beta^8 = N^2*alpha^2, N = w^2
 */
function [3:0] GF_SQ_SCL_4;
  input [3:0] A;
  reg [1:0] a, b, ab2, b2, b2N2;
  begin
    a = A[3:2];
    b = A[1:0];
    ab2 = GF_SQ_2(a ^ b);
    b2 = GF_SQ_2(b);
    b2N2 = GF_SCLW_2(b2);

    GF_SQ_SCL_4 = { ab2, b2N2 };
  end
endfunction

/* multiply in GF(2^4)/GF(2^2), shared factors, basis [alpha^8, alpha^2] */
function [3:0] GF_MULS_4;
  input [3:0] A;
  input [1:0] a;
  input Al;
  input Ah;
  input aa;
  input [3:0] B;
  input [1:0] b;
  input Bl;
  input Bh;
  input bb;
  reg [1:0] ph, pl, ps, p;
  reg t;
  begin

    ph = GF_MULS_2 (A[3:2], Ah, B[3:2], Bh);
    pl = GF_MULS_2 (A[1:0], Al, B[1:0], Bl);
    p = GF_MULS_SCL_2 ( a, aa, b, bb);
    GF_MULS_4 = { (ph ^ p), (pl ^ p) };
  end
endfunction

/* inverse in GF(2^8)/GF(2^4), using normal basis [d^16, d] */
function [7:0] GF_INV_8;
  input [7:0] A;
  reg [3:0] a, b, c, d, p, q;
  reg [1:0] sa, sb, sd, t; /* for shared factors in multipliers */ 
  reg al, ah, aa, bl, bh, bb, dl, dh, dd; /* for shared factors */
  reg c1, c2, c3; /* for temp var */
  begin

  a = A[7:4];
  b = A[3:0];
  sa = a[3:2] ^ a[1:0];
  sb = b[3:2] ^ b[1:0];
  al = a[1] ^ a[0];
  ah = a[3] ^ a[2];
  aa = sa[1] ^ sa[0];
  bl = b[1] ^ b[0];
  bh = b[3] ^ b[2];
  bb = sb[1] ^ sb[0];
  
  /* optimize this section as shown below
  GF_MULS_4 abmul(a, sa, al, ah, aa, b, sb, bl, bh, bb, ab);
  GF_SQ_SCL_4 absq( (a ^ b), ab2);
  GF_INV_4 dinv( (ab ^ ab2), d);
  */
    c1 = ~(ah & bh);
    c2 = ~(sa[0] & sb[0]);
    c3 = ~(aa & bb);
    c = {
    (~(sa[0] | sb[0]) ^ (~(a[3] & b[3]))) ^ c1 ^ c3 ,
    (~(sa[1] | sb[1]) ^ (~(a[2] & b[2]))) ^ c1 ^ c2 ,
    (~(al | bl) ^ (~(a[1] & b[1]))) ^ c2 ^ c3 ,
    (~(a[0] | b[0]) ^ (~(al & bl))) ^ (~(sa[1] & sb[1])) ^ c2
    };
    d = GF_INV_4( c);
  /* end of optimization */
  sd = d[3:2] ^ d[1:0];
  dl = d[1] ^ d[0];
  dh = d[3] ^ d[2];
  dd = sd[1] ^ sd[0];
  p = GF_MULS_4(d, sd, dl, dh, dd, b, sb, bl, bh, bb);
  q = GF_MULS_4(d, sd, dl, dh, dd, a, sa, al, ah, aa);
  GF_INV_8 = { p, q };
  end
endfunction

/* select and invert (NOT) byte, using MUX21I */
function [7:0] SELECT_NOT_8;
  input [7:0] A;
  input [7:0] B;
  input s;
  begin
    SELECT_NOT_8 = ~( s ? A : B );
  end
endfunction


  wire [7:0] B, C, D, X, Y, Z;
  wire R1, R2, R3, R4, R5, R6, R7, R8, R9;
  wire T1, T2, T3, T4, T5, T6, T7, T8, T9, T10;

  /* change basis from GF(2^8)
  /* combine with bit inverse matrix multiply of Sbox */
  assign R1 = A[7] ^ A[5] ;
  assign R2 = A[7] ~^ A[4] ;
  assign R3 = A[6] ^ A[0] ;
  assign R4 = A[5] ~^ R3 ;
  assign R5 = A[4] ^ R4 ;
  assign R6 = A[3] ^ A[0] ;
  assign R7 = A[2] ^ R1 ;
  assign R8 = A[1] ^ R3 ;
  assign R9 = A[3] ^ R8 ;
  assign B[7] = R7 ~^ R8 ;
  assign B[6] = R5 ;
  assign B[5] = A[1] ^ R4 ;
  assign B[4] = R1 ~^ R3 ;
  assign B[3] = A[1]^ R2 ^ R6 ;
  assign B[2] = ~ A[0] ;
  assign B[1] = R4 ;
  assign B[0] = A[2] ~^ R9 ;
  assign Y[7] = R2 ;
  assign Y[6] = A[4] ^ R8 ;
  
  assign Y[5] = A[6] ^ A[4] ;
  assign Y[4] = R9 ;
  assign Y[3] = A[6] ~^ R2 ;
  assign Y[2] = R7 ;
  assign Y[1] = A[4] ^ R6 ;
  assign Y[0] = A[1] ^ R5 ;
  
  assign Z = SELECT_NOT_8 ( B, Y, encrypt );
  assign C = GF_INV_8( Z );
  /* change basis back from GF(2^8)/GF(2^4)/GF(2^2) to GF(2^8) */
  assign T1 = C[7] ^ C[3];
  assign T2 = C[6] ^ C[4];
  assign T3 = C[6] ^ C[0];
  assign T4 = C[5] ~^ C[3] ; 
  assign T5 = C[5] ~^ T1 ; 
  assign T6 = C[5] ~^ C[1] ; 
  assign T7 = C[4] ~^ T6 ; 
  assign T8 = C[2] ^ T4 ; 
  assign T9 = C[1] ^ T2 ; 
  assign T10 = T3 ^ T5 ; 
  assign D[7] = T4 ;
  assign D[6] = T1 ;
  assign D[5] = T3 ;
  assign D[4] = T5 ;
  assign D[3] = T2 ^ T5;
  assign D[2] = T3 ^ T8;
  assign D[1] = T7 ;
  assign D[0] = T9 ;
  assign X[7] = C[4] ~^ C[1] ;
  assign X[6] = C[1] ^ T10 ;
  assign X[5] = C[2] ^ T10 ;
  assign X[4] = C[6] ~^ C[1] ;
  assign X[3] = T8 ^ T9 ;
  assign X[2] = C[7] ~^ T7 ;
  assign X[1] = T6 ;
  assign X[0] = ~ C[2] ;
  assign Q = SELECT_NOT_8 ( D, X, encrypt );
endmodule
