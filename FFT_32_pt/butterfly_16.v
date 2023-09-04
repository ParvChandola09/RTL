module butterfly_16 #
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
	output signed  [N-1:0]  y15_i   
	);


//1st butterfly
butterfly_2 # (.N(N), .Q(Q)) B1
	(
	.x0_r(x0_r),
	.x0_i(x0_i),
	.x1_r(x8_r),
	.x1_i(x8_i),
	.tw_r(tw1_r),
	.tw_i(tw1_i),
	.y0_r(y0_r),
	.y0_i(y0_i),
	.y1_r(y8_r),
	.y1_i(y8_i)
	);


//2nd butterfly
butterfly_2 # (.N(N), .Q(Q)) B2
	(
	.x0_r(x1_r),
	.x0_i(x1_i),
	.x1_r(x9_r),
	.x1_i(x9_i),
	.tw_r(tw2_r),
	.tw_i(tw2_i),
	.y0_r(y1_r),
	.y0_i(y1_i),
	.y1_r(y9_r),
	.y1_i(y9_i)
	);

//3rd butterfly
butterfly_2 # (.N(N), .Q(Q)) B3
	(
	.x0_r(x2_r),
	.x0_i(x2_i),
	.x1_r(x10_r),
	.x1_i(x10_i),
	.tw_r(tw3_r),
	.tw_i(tw3_i),
	.y0_r(y2_r),
	.y0_i(y2_i),
	.y1_r(y10_r),
	.y1_i(y10_i)
	);

//4th butterfly
butterfly_2 # (.N(N), .Q(Q)) B4
	(
	.x0_r(x3_r),
	.x0_i(x3_i),
	.x1_r(x11_r),
	.x1_i(x11_i),
	.tw_r(tw4_r),
	.tw_i(tw4_i),
	.y0_r(y3_r),
	.y0_i(y3_i),
	.y1_r(y11_r),
	.y1_i(y11_i)
	);

//5th butterfly
butterfly_2 # (.N(N), .Q(Q)) B5
	(
	.x0_r(x4_r),
	.x0_i(x4_i),
	.x1_r(x12_r),
	.x1_i(x12_i),
	.tw_r(tw5_r),
	.tw_i(tw5_i),
	.y0_r(y4_r),
	.y0_i(y4_i),
	.y1_r(y12_r),
	.y1_i(y12_i)
	);

//6th butterfly
butterfly_2 # (.N(N), .Q(Q)) B6
	(
	.x0_r(x5_r),
	.x0_i(x5_i),
	.x1_r(x13_r),
	.x1_i(x13_i),
	.tw_r(tw6_r),
	.tw_i(tw6_i),
	.y0_r(y5_r),
	.y0_i(y5_i),
	.y1_r(y13_r),
	.y1_i(y13_i)
	);

//7th butterfly
butterfly_2 # (.N(N), .Q(Q)) B7
	(
	.x0_r(x6_r),
	.x0_i(x6_i),
	.x1_r(x14_r),
	.x1_i(x14_i),
	.tw_r(tw7_r),
	.tw_i(tw7_i),
	.y0_r(y6_r),
	.y0_i(y6_i),
	.y1_r(y14_r),
	.y1_i(y14_i)
	);

//8th butterfly
butterfly_2 # (.N(N), .Q(Q)) B8
	(
	.x0_r(x7_r),
	.x0_i(x7_i),
	.x1_r(x15_r),
	.x1_i(x15_i),
	.tw_r(tw8_r),
	.tw_i(tw8_i),
	.y0_r(y7_r),
	.y0_i(y7_i),
	.y1_r(y15_r),
	.y1_i(y15_i)
	);

endmodule