module butterfly_4 #
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

	input  signed  [N-1:0]  tw1_r,
	input  signed  [N-1:0]  tw1_i,
	input  signed  [N-1:0]  tw2_r,
	input  signed  [N-1:0]  tw2_i,

	output signed  [N-1:0]  y0_r, 
	output signed  [N-1:0]  y0_i,
	output signed  [N-1:0]  y1_r,
	output signed  [N-1:0]  y1_i,
	output signed  [N-1:0]  y2_r,
	output signed  [N-1:0]  y2_i, 
	output signed  [N-1:0]  y3_r,
	output signed  [N-1:0]  y3_i   
	);


//1st butterfly
butterfly_2 # (.N(N), .Q(Q)) B1
	(
	.x0_r(x0_r),
	.x0_i(x0_i),
	.x1_r(x2_r),
	.x1_i(x2_i),
	.tw_r(tw1_r),
	.tw_i(tw1_i),
	.y0_r(y0_r),
	.y0_i(y0_i),
	.y1_r(y2_r),
	.y1_i(y2_i)
	);


//2nd butterfly
butterfly_2 # (.N(N), .Q(Q)) B2
	(
	.x0_r(x1_r),
	.x0_i(x1_i),
	.x1_r(x3_r),
	.x1_i(x3_i),
	.tw_r(tw2_r),
	.tw_i(tw2_i),
	.y0_r(y1_r),
	.y0_i(y1_i),
	.y1_r(y3_r),
	.y1_i(y3_i)
	);


endmodule
