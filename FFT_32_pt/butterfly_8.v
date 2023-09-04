module butterfly_8 #
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

	input  signed  [N-1:0]  tw1_r,
	input  signed  [N-1:0]  tw1_i,
	input  signed  [N-1:0]  tw2_r,
	input  signed  [N-1:0]  tw2_i,
	input  signed  [N-1:0]  tw3_r,
	input  signed  [N-1:0]  tw3_i,
	input  signed  [N-1:0]  tw4_r,
	input  signed  [N-1:0]  tw4_i,

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
	output signed  [N-1:0]  y7_i   
	);


//1st butterfly
butterfly_2 # (.N(N), .Q(Q)) B1
	(
	.x0_r(x0_r),
	.x0_i(x0_i),
	.x1_r(x4_r),
	.x1_i(x4_i),
	.tw_r(tw1_r),
	.tw_i(tw1_i),
	.y0_r(y0_r),
	.y0_i(y0_i),
	.y1_r(y4_r),
	.y1_i(y4_i)
	);


//2nd butterfly
butterfly_2 # (.N(N), .Q(Q)) B2
	(
	.x0_r(x1_r),
	.x0_i(x1_i),
	.x1_r(x5_r),
	.x1_i(x5_i),
	.tw_r(tw2_r),
	.tw_i(tw2_i),
	.y0_r(y1_r),
	.y0_i(y1_i),
	.y1_r(y5_r),
	.y1_i(y5_i)
	);

//3rd butterfly
butterfly_2 # (.N(N), .Q(Q)) B3
	(
	.x0_r(x2_r),
	.x0_i(x2_i),
	.x1_r(x6_r),
	.x1_i(x6_i),
	.tw_r(tw3_r),
	.tw_i(tw3_i),
	.y0_r(y2_r),
	.y0_i(y2_i),
	.y1_r(y6_r),
	.y1_i(y6_i)
	);

//4th butterfly
butterfly_2 # (.N(N), .Q(Q)) B4
	(
	.x0_r(x3_r),
	.x0_i(x3_i),
	.x1_r(x7_r),
	.x1_i(x7_i),
	.tw_r(tw4_r),
	.tw_i(tw4_i),
	.y0_r(y3_r),
	.y0_i(y3_i),
	.y1_r(y7_r),
	.y1_i(y7_i)
	);

endmodule
