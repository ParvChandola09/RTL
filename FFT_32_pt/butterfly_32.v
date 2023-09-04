module butterfly_32 #
	(
	parameter N = 16,
	parameter Q = 8
	)
	(
	input  signed  [N-1:0]  x0_r,
	input  signed  [N-1:0]  x0_i,
	input  signed  [N-1:0]  x1_r,
	input  signed  [N-1:0]  x1_i,
	input  signed  [N-1:0]  x2_r,
	input  signed  [N-1:0]  x2_i,
	input  signed  [N-1:0]  x3_r,
	input  signed  [N-1:0]  x3_i,
	input  signed  [N-1:0]  x4_r,
	input  signed  [N-1:0]  x4_i,
	input  signed  [N-1:0]  x5_r,
	input  signed  [N-1:0]  x5_i,
	input  signed  [N-1:0]  x6_r,
	input  signed  [N-1:0]  x6_i,
	input  signed  [N-1:0]  x7_r,
	input  signed  [N-1:0]  x7_i,
	input  signed  [N-1:0]  x8_r,
	input  signed  [N-1:0]  x8_i,
	input  signed  [N-1:0]  x9_r,
	input  signed  [N-1:0]  x9_i,
	input  signed  [N-1:0]  x10_r,
	input  signed  [N-1:0]  x10_i,
	input  signed  [N-1:0]  x11_r,
	input  signed  [N-1:0]  x11_i,
	input  signed  [N-1:0]  x12_r,
	input  signed  [N-1:0]  x12_i,
	input  signed  [N-1:0]  x13_r,
	input  signed  [N-1:0]  x13_i,
	input  signed  [N-1:0]  x14_r,
	input  signed  [N-1:0]  x14_i,
	input  signed  [N-1:0]  x15_r,
	input  signed  [N-1:0]  x15_i,
	input  signed  [N-1:0]  x16_r,
	input  signed  [N-1:0]  x16_i,
	input  signed  [N-1:0]  x17_r,
	input  signed  [N-1:0]  x17_i,
	input  signed  [N-1:0]  x18_r,
	input  signed  [N-1:0]  x18_i,
	input  signed  [N-1:0]  x19_r,
	input  signed  [N-1:0]  x19_i,
	input  signed  [N-1:0]  x20_r,
	input  signed  [N-1:0]  x20_i,
	input  signed  [N-1:0]  x21_r,
	input  signed  [N-1:0]  x21_i,
	input  signed  [N-1:0]  x22_r,
	input  signed  [N-1:0]  x22_i,
	input  signed  [N-1:0]  x23_r,
	input  signed  [N-1:0]  x23_i,
	input  signed  [N-1:0]  x24_r,
	input  signed  [N-1:0]  x24_i,
	input  signed  [N-1:0]  x25_r,
	input  signed  [N-1:0]  x25_i,
	input  signed  [N-1:0]  x26_r,
	input  signed  [N-1:0]  x26_i,
	input  signed  [N-1:0]  x27_r,
	input  signed  [N-1:0]  x27_i,
	input  signed  [N-1:0]  x28_r,
	input  signed  [N-1:0]  x28_i,
	input  signed  [N-1:0]  x29_r,
	input  signed  [N-1:0]  x29_i,
	input  signed  [N-1:0]  x30_r,
	input  signed  [N-1:0]  x30_i,
	input  signed  [N-1:0]  x31_r,
	input  signed  [N-1:0]  x31_i,

	input  signed  [N-1:0]  tw1_r,
	input  signed  [N-1:0]  tw1_i,
	input  signed  [N-1:0]  tw2_r,
	input  signed  [N-1:0]  tw2_i,
	input  signed  [N-1:0]  tw3_r,
	input  signed  [N-1:0]  tw3_i,
	input  signed  [N-1:0]  tw4_r,
	input  signed  [N-1:0]  tw4_i,
	input  signed  [N-1:0]  tw5_r,
	input  signed  [N-1:0]  tw5_i,
	input  signed  [N-1:0]  tw6_r,
	input  signed  [N-1:0]  tw6_i,
	input  signed  [N-1:0]  tw7_r,
	input  signed  [N-1:0]  tw7_i,
	input  signed  [N-1:0]  tw8_r,
	input  signed  [N-1:0]  tw8_i,
	input  signed  [N-1:0]  tw9_r,
	input  signed  [N-1:0]  tw9_i,
	input  signed  [N-1:0]  tw10_r,
	input  signed  [N-1:0]  tw10_i,
	input  signed  [N-1:0]  tw11_r,
	input  signed  [N-1:0]  tw11_i,
	input  signed  [N-1:0]  tw12_r,
	input  signed  [N-1:0]  tw12_i,
	input  signed  [N-1:0]  tw13_r,
	input  signed  [N-1:0]  tw13_i,
	input  signed  [N-1:0]  tw14_r,
	input  signed  [N-1:0]  tw14_i,
	input  signed  [N-1:0]  tw15_r,
	input  signed  [N-1:0]  tw15_i,
	input  signed  [N-1:0]  tw16_r,
	input  signed  [N-1:0]  tw16_i,

	output signed  [N-1:0]  y0_r, 
	output signed  [N-1:0]  y0_i,
	output signed  [N-1:0]  y1_r,
	output signed  [N-1:0]  y1_i,
	output signed  [N-1:0]  y2_r,
	output signed  [N-1:0]  y2_i, 
	output signed  [N-1:0]  y3_r,
	output signed  [N-1:0]  y3_i, 
	output signed  [N-1:0]  y4_r, 
	output signed  [N-1:0]  y4_i,
	output signed  [N-1:0]  y5_r,
	output signed  [N-1:0]  y5_i,
	output signed  [N-1:0]  y6_r,
	output signed  [N-1:0]  y6_i, 
	output signed  [N-1:0]  y7_r,
	output signed  [N-1:0]  y7_i,  
	output signed  [N-1:0]  y8_r, 
	output signed  [N-1:0]  y8_i,
	output signed  [N-1:0]  y9_r,
	output signed  [N-1:0]  y9_i,
	output signed  [N-1:0]  y10_r,
	output signed  [N-1:0]  y10_i, 
	output signed  [N-1:0]  y11_r,
	output signed  [N-1:0]  y11_i, 
	output signed  [N-1:0]  y12_r, 
	output signed  [N-1:0]  y12_i,
	output signed  [N-1:0]  y13_r,
	output signed  [N-1:0]  y13_i,
	output signed  [N-1:0]  y14_r,
	output signed  [N-1:0]  y14_i, 
	output signed  [N-1:0]  y15_r,
	output signed  [N-1:0]  y15_i,  
	output signed  [N-1:0]  y16_r, 
	output signed  [N-1:0]  y16_i,
	output signed  [N-1:0]  y17_r,
	output signed  [N-1:0]  y17_i,
	output signed  [N-1:0]  y18_r,
	output signed  [N-1:0]  y18_i, 
	output signed  [N-1:0]  y19_r,
	output signed  [N-1:0]  y19_i, 
	output signed  [N-1:0]  y20_r, 
	output signed  [N-1:0]  y20_i,
	output signed  [N-1:0]  y21_r,
	output signed  [N-1:0]  y21_i,
	output signed  [N-1:0]  y22_r,
	output signed  [N-1:0]  y22_i, 
	output signed  [N-1:0]  y23_r,
	output signed  [N-1:0]  y23_i,  
	output signed  [N-1:0]  y24_r, 
	output signed  [N-1:0]  y24_i,
	output signed  [N-1:0]  y25_r,
	output signed  [N-1:0]  y25_i,
	output signed  [N-1:0]  y26_r,
	output signed  [N-1:0]  y26_i, 
	output signed  [N-1:0]  y27_r,
	output signed  [N-1:0]  y27_i, 
	output signed  [N-1:0]  y28_r, 
	output signed  [N-1:0]  y28_i,
	output signed  [N-1:0]  y29_r,
	output signed  [N-1:0]  y29_i,
	output signed  [N-1:0]  y30_r,
	output signed  [N-1:0]  y30_i, 
	output signed  [N-1:0]  y31_r,
	output signed  [N-1:0]  y31_i   
	);


//1st butterfly
butterfly_2 # (.N(N), .Q(Q)) B1
	(
	.x0_r(x0_r),
	.x0_i(x0_i),
	.x1_r(x16_r),
	.x1_i(x16_i),
	.tw_r(tw1_r),
	.tw_i(tw1_i),
	.y0_r(y0_r),
	.y0_i(y0_i),
	.y1_r(y16_r),
	.y1_i(y16_i)
	);


//2nd butterfly
butterfly_2 # (.N(N), .Q(Q)) B2
	(
	.x0_r(x1_r),
	.x0_i(x1_i),
	.x1_r(x17_r),
	.x1_i(x17_i),
	.tw_r(tw2_r),
	.tw_i(tw2_i),
	.y0_r(y1_r),
	.y0_i(y1_i),
	.y1_r(y17_r),
	.y1_i(y17_i)
	);

//3rd butterfly
butterfly_2 # (.N(N), .Q(Q)) B3
	(
	.x0_r(x2_r),
	.x0_i(x2_i),
	.x1_r(x18_r),
	.x1_i(x18_i),
	.tw_r(tw3_r),
	.tw_i(tw3_i),
	.y0_r(y2_r),
	.y0_i(y2_i),
	.y1_r(y18_r),
	.y1_i(y18_i)
	);

//4th butterfly
butterfly_2 # (.N(N), .Q(Q)) B4
	(
	.x0_r(x3_r),
	.x0_i(x3_i),
	.x1_r(x19_r),
	.x1_i(x19_i),
	.tw_r(tw4_r),
	.tw_i(tw4_i),
	.y0_r(y3_r),
	.y0_i(y3_i),
	.y1_r(y19_r),
	.y1_i(y19_i)
	);

//5th butterfly
butterfly_2 # (.N(N), .Q(Q)) B5
	(
	.x0_r(x4_r),
	.x0_i(x4_i),
	.x1_r(x20_r),
	.x1_i(x20_i),
	.tw_r(tw5_r),
	.tw_i(tw5_i),
	.y0_r(y4_r),
	.y0_i(y4_i),
	.y1_r(y20_r),
	.y1_i(y20_i)
	);

//6th butterfly
butterfly_2 # (.N(N), .Q(Q)) B6
	(
	.x0_r(x5_r),
	.x0_i(x5_i),
	.x1_r(x21_r),
	.x1_i(x21_i),
	.tw_r(tw6_r),
	.tw_i(tw6_i),
	.y0_r(y5_r),
	.y0_i(y5_i),
	.y1_r(y21_r),
	.y1_i(y21_i)
	);

//7th butterfly
butterfly_2 # (.N(N), .Q(Q)) B7
	(
	.x0_r(x6_r),
	.x0_i(x6_i),
	.x1_r(x22_r),
	.x1_i(x22_i),
	.tw_r(tw7_r),
	.tw_i(tw7_i),
	.y0_r(y6_r),
	.y0_i(y6_i),
	.y1_r(y22_r),
	.y1_i(y22_i)
	);

//8th butterfly
butterfly_2 # (.N(N), .Q(Q)) B8
	(
	.x0_r(x7_r),
	.x0_i(x7_i),
	.x1_r(x23_r),
	.x1_i(x23_i),
	.tw_r(tw8_r),
	.tw_i(tw8_i),
	.y0_r(y7_r),
	.y0_i(y7_i),
	.y1_r(y23_r),
	.y1_i(y23_i)
	);

//9th butterfly
butterfly_2 # (.N(N), .Q(Q)) B9
	(
	.x0_r(x8_r),
	.x0_i(x8_i),
	.x1_r(x24_r),
	.x1_i(x24_i),
	.tw_r(tw9_r),
	.tw_i(tw9_i),
	.y0_r(y8_r),
	.y0_i(y8_i),
	.y1_r(y24_r),
	.y1_i(y24_i)
	);

//10th butterfly
butterfly_2 # (.N(N), .Q(Q)) B10
	(
	.x0_r(x9_r),
	.x0_i(x9_i),
	.x1_r(x25_r),
	.x1_i(x25_i),
	.tw_r(tw10_r),
	.tw_i(tw10_i),
	.y0_r(y9_r),
	.y0_i(y9_i),
	.y1_r(y25_r),
	.y1_i(y25_i)
	);

//11th butterfly
butterfly_2 # (.N(N), .Q(Q)) B11
	(
	.x0_r(x10_r),
	.x0_i(x10_i),
	.x1_r(x26_r),
	.x1_i(x26_i),
	.tw_r(tw11_r),
	.tw_i(tw11_i),
	.y0_r(y10_r),
	.y0_i(y10_i),
	.y1_r(y26_r),
	.y1_i(y26_i)
	);

//12th butterfly
butterfly_2 # (.N(N), .Q(Q)) B12
	(
	.x0_r(x11_r),
	.x0_i(x11_i),
	.x1_r(x27_r),
	.x1_i(x27_i),
	.tw_r(tw12_r),
	.tw_i(tw12_i),
	.y0_r(y11_r),
	.y0_i(y11_i),
	.y1_r(y27_r),
	.y1_i(y27_i)
	);

//13th butterfly
butterfly_2 # (.N(N), .Q(Q)) B13
	(
	.x0_r(x12_r),
	.x0_i(x12_i),
	.x1_r(x28_r),
	.x1_i(x28_i),
	.tw_r(tw13_r),
	.tw_i(tw13_i),
	.y0_r(y12_r),
	.y0_i(y12_i),
	.y1_r(y28_r),
	.y1_i(y28_i)
	);

//14th butterfly
butterfly_2 # (.N(N), .Q(Q)) B14
	(
	.x0_r(x13_r),
	.x0_i(x13_i),
	.x1_r(x29_r),
	.x1_i(x29_i),
	.tw_r(tw14_r),
	.tw_i(tw14_i),
	.y0_r(y13_r),
	.y0_i(y13_i),
	.y1_r(y29_r),
	.y1_i(y29_i)
	);

//15th butterfly
butterfly_2 # (.N(N), .Q(Q)) B15
	(
	.x0_r(x14_r),
	.x0_i(x14_i),
	.x1_r(x30_r),
	.x1_i(x30_i),
	.tw_r(tw15_r),
	.tw_i(tw15_i),
	.y0_r(y14_r),
	.y0_i(y14_i),
	.y1_r(y30_r),
	.y1_i(y30_i)
	);

//16th butterfly
butterfly_2 # (.N(N), .Q(Q)) B16
	(
	.x0_r(x15_r),
	.x0_i(x15_i),
	.x1_r(x31_r),
	.x1_i(x31_i),
	.tw_r(tw16_r),
	.tw_i(tw16_i),
	.y0_r(y15_r),
	.y0_i(y15_i),
	.y1_r(y31_r),
	.y1_i(y31_i)
	);

endmodule
