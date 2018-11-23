module reference_sbox
  (
   input [7:0]      index,
   output reg [7:0] sbox
   );
   
  always @*
    begin
      case (index)
	0 : sbox = 8'h63;
	1 : sbox = 8'h7c;
	2 : sbox = 8'h77;
	3 : sbox = 8'h7b;
	4 : sbox = 8'hf2;
	5 : sbox = 8'h6b;
	6 : sbox = 8'h6f;
	7 : sbox = 8'hc5;
	8 : sbox = 8'h30;
	9 : sbox = 8'h01;
	10 : sbox = 8'h67;
	11 : sbox = 8'h2b;
	12 : sbox = 8'hfe;
	13 : sbox = 8'hd7;
	14 : sbox = 8'hab;
	15 : sbox = 8'h76;
	16 : sbox = 8'hca;
	17 : sbox = 8'h82;
	18 : sbox = 8'hc9;
	19 : sbox = 8'h7d;
	20 : sbox = 8'hfa;
	21 : sbox = 8'h59;
	22 : sbox = 8'h47;
	23 : sbox = 8'hf0;
	24 : sbox = 8'had;
	25 : sbox = 8'hd4;
	26 : sbox = 8'ha2;
	27 : sbox = 8'haf;
	28 : sbox = 8'h9c;
	29 : sbox = 8'ha4;
	30 : sbox = 8'h72;
	31 : sbox = 8'hc0;
	32 : sbox = 8'hb7;
	33 : sbox = 8'hfd;
	34 : sbox = 8'h93;
	35 : sbox = 8'h26;
	36 : sbox = 8'h36;
	37 : sbox = 8'h3f;
	38 : sbox = 8'hf7;
	39 : sbox = 8'hcc;
	40 : sbox = 8'h34;
	41 : sbox = 8'ha5;
	42 : sbox = 8'he5;
	43 : sbox = 8'hf1;
	44 : sbox = 8'h71;
	45 : sbox = 8'hd8;
	46 : sbox = 8'h31;
	47 : sbox = 8'h15;
	48 : sbox = 8'h04;
	49 : sbox = 8'hc7;
	50 : sbox = 8'h23;
	51 : sbox = 8'hc3;
	52 : sbox = 8'h18;
	53 : sbox = 8'h96;
	54 : sbox = 8'h05;
	55 : sbox = 8'h9a;
	56 : sbox = 8'h07;
	57 : sbox = 8'h12;
	58 : sbox = 8'h80;
	59 : sbox = 8'he2;
	60 : sbox = 8'heb;
	61 : sbox = 8'h27;
	62 : sbox = 8'hb2;
	63 : sbox = 8'h75;
	64 : sbox = 8'h09;
	65 : sbox = 8'h83;
	66 : sbox = 8'h2c;
	67 : sbox = 8'h1a;
	68 : sbox = 8'h1b;
	69 : sbox = 8'h6e;
	70 : sbox = 8'h5a;
	71 : sbox = 8'ha0;
	72 : sbox = 8'h52;
	73 : sbox = 8'h3b;
	74 : sbox = 8'hd6;
	75 : sbox = 8'hb3;
	76 : sbox = 8'h29;
	77 : sbox = 8'he3;
	78 : sbox = 8'h2f;
	79 : sbox = 8'h84;
	80 : sbox = 8'h53;
	81 : sbox = 8'hd1;
	82 : sbox = 8'h00;
	83 : sbox = 8'hed;
	84 : sbox = 8'h20;
	85 : sbox = 8'hfc;
	86 : sbox = 8'hb1;
	87 : sbox = 8'h5b;
	88 : sbox = 8'h6a;
	89 : sbox = 8'hcb;
	90 : sbox = 8'hbe;
	91 : sbox = 8'h39;
	92 : sbox = 8'h4a;
	93 : sbox = 8'h4c;
	94 : sbox = 8'h58;
	95 : sbox = 8'hcf;
	96 : sbox = 8'hd0;
	97 : sbox = 8'hef;
	98 : sbox = 8'haa;
	99 : sbox = 8'hfb;
	100 : sbox = 8'h43;
	101 : sbox = 8'h4d;
	102 : sbox = 8'h33;
	103 : sbox = 8'h85;
	104 : sbox = 8'h45;
	105 : sbox = 8'hf9;
	106 : sbox = 8'h02;
	107 : sbox = 8'h7f;
	108 : sbox = 8'h50;
	109 : sbox = 8'h3c;
	110 : sbox = 8'h9f;
	111 : sbox = 8'ha8;
	112 : sbox = 8'h51;
	113 : sbox = 8'ha3;
	114 : sbox = 8'h40;
	115 : sbox = 8'h8f;
	116 : sbox = 8'h92;
	117 : sbox = 8'h9d;
	118 : sbox = 8'h38;
	119 : sbox = 8'hf5;
	120 : sbox = 8'hbc;
	121 : sbox = 8'hb6;
	122 : sbox = 8'hda;
	123 : sbox = 8'h21;
	124 : sbox = 8'h10;
	125 : sbox = 8'hff;
	126 : sbox = 8'hf3;
	127 : sbox = 8'hd2;
	128 : sbox = 8'hcd;
	129 : sbox = 8'h0c;
	130 : sbox = 8'h13;
	131 : sbox = 8'hec;
	132 : sbox = 8'h5f;
	133 : sbox = 8'h97;
	134 : sbox = 8'h44;
	135 : sbox = 8'h17;
	136 : sbox = 8'hc4;
	137 : sbox = 8'ha7;
	138 : sbox = 8'h7e;
	139 : sbox = 8'h3d;
	140 : sbox = 8'h64;
	141 : sbox = 8'h5d;
	142 : sbox = 8'h19;
	143 : sbox = 8'h73;
	144 : sbox = 8'h60;
	145 : sbox = 8'h81;
	146 : sbox = 8'h4f;
	147 : sbox = 8'hdc;
	148 : sbox = 8'h22;
	149 : sbox = 8'h2a;
	150 : sbox = 8'h90;
	151 : sbox = 8'h88;
	152 : sbox = 8'h46;
	153 : sbox = 8'hee;
	154 : sbox = 8'hb8;
	155 : sbox = 8'h14;
	156 : sbox = 8'hde;
	157 : sbox = 8'h5e;
	158 : sbox = 8'h0b;
	159 : sbox = 8'hdb;
	160 : sbox = 8'he0;
	161 : sbox = 8'h32;
	162 : sbox = 8'h3a;
	163 : sbox = 8'h0a;
	164 : sbox = 8'h49;
	165 : sbox = 8'h06;
	166 : sbox = 8'h24;
	167 : sbox = 8'h5c;
	168 : sbox = 8'hc2;
	169 : sbox = 8'hd3;
	170 : sbox = 8'hac;
	171 : sbox = 8'h62;
	172 : sbox = 8'h91;
	173 : sbox = 8'h95;
	174 : sbox = 8'he4;
	175 : sbox = 8'h79;
	176 : sbox = 8'he7;
	177 : sbox = 8'hc8;
	178 : sbox = 8'h37;
	179 : sbox = 8'h6d;
	180 : sbox = 8'h8d;
	181 : sbox = 8'hd5;
	182 : sbox = 8'h4e;
	183 : sbox = 8'ha9;
	184 : sbox = 8'h6c;
	185 : sbox = 8'h56;
	186 : sbox = 8'hf4;
	187 : sbox = 8'hea;
	188 : sbox = 8'h65;
	189 : sbox = 8'h7a;
	190 : sbox = 8'hae;
	191 : sbox = 8'h8;
	192 : sbox = 8'hba;
	193 : sbox = 8'h78;
	194 : sbox = 8'h25;
	195 : sbox = 8'h2e;
	196 : sbox = 8'h1c;
	197 : sbox = 8'ha6;
	198 : sbox = 8'hb4;
	199 : sbox = 8'hc6;
	200 : sbox = 8'he8;
	201 : sbox = 8'hdd;
	202 : sbox = 8'h74;
	203 : sbox = 8'h1f;
	204 : sbox = 8'h4b;
	205 : sbox = 8'hbd;
	206 : sbox = 8'h8b;
	207 : sbox = 8'h8a;
	208 : sbox = 8'h70;
	209 : sbox = 8'h3e;
	210 : sbox = 8'hb5;
	211 : sbox = 8'h66;
	212 : sbox = 8'h48;
	213 : sbox = 8'h03;
	214 : sbox = 8'hf6;
	215 : sbox = 8'he;
	216 : sbox = 8'h61;
	217 : sbox = 8'h35;
	218 : sbox = 8'h57;
	219 : sbox = 8'hb9;
	220 : sbox = 8'h86;
	221 : sbox = 8'hc1;
	222 : sbox = 8'h1d;
	223 : sbox = 8'h9e;
	224 : sbox = 8'he1;
	225 : sbox = 8'hf8;
	226 : sbox = 8'h98;
	227 : sbox = 8'h11;
	228 : sbox = 8'h69;
	229 : sbox = 8'hd9;
	230 : sbox = 8'h8e;
	231 : sbox = 8'h94;
	232 : sbox = 8'h9b;
	233 : sbox = 8'h1e;
	234 : sbox = 8'h87;
	235 : sbox = 8'he9;
	236 : sbox = 8'hce;
	237 : sbox = 8'h55;
	238 : sbox = 8'h28;
	239 : sbox = 8'hdf;
	240 : sbox = 8'h8c;
	241 : sbox = 8'ha1;
	242 : sbox = 8'h89;
	243 : sbox = 8'h0d;
	244 : sbox = 8'hbf;
	245 : sbox = 8'he6;
	246 : sbox = 8'h42;
	247 : sbox = 8'h68;
	248 : sbox = 8'h41;
	249 : sbox = 8'h99;
	250 : sbox = 8'h2d;
	251 : sbox = 8'h0f;
	252 : sbox = 8'hb0;
	253 : sbox = 8'h54;
	254 : sbox = 8'hbb;
	255 : sbox = 8'h16;
      endcase // case(index)
    end // always @ *

endmodule 
