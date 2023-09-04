module butterfly_2 #
	(
	parameter N = 16,
	parameter Q = 8
	)
	(
	input  signed  [N-1:0]  x0_r,
	input  signed  [N-1:0]  x0_i,
	input  signed  [N-1:0]  x1_r,
	input  signed  [N-1:0]  x1_i,

	input  signed  [N-1:0]  tw_r,
	input  signed  [N-1:0]  tw_i,

	output  signed  [N-1:0]  y0_r,
	output 	signed  [N-1:0]  y0_i,
	output  signed  [N-1:0]  y1_r,
	output  signed  [N-1:0]  y1_i
	);

//calculation of the multiplication stage(multiplication with twiddle factors)

wire  signed  [N-1:0]  x1_mul_R;
wire  signed  [N-1:0]  x1_mul_I;


//multiplication with the twiddle factor
complex_multiplier # (.N(N), .Q(Q)) M1
	(
	.x1_real(x1_r),
	.x1_imag(x1_i),
	.x2_real(tw_r),
	.x2_imag(tw_i),
	.z_real(x1_mul_R),
	.z_imag(x1_mul_I)
	);   

//1st equation output 
complex_adder # (.N(N)) A1
	(
	.x_real(x0_r),
	.x_imag(x0_i),
	.y_real(x1_mul_R),
	.y_imag(x1_mul_I),
	.z_real(y0_r),
	.z_imag(y0_i)
	);

//2nd equation output
complex_subtractor # (.N(N)) S1
	(
	.x_real(x0_r),
	.x_imag(x0_i),
	.y_real(x1_mul_R),
	.y_imag(x1_mul_I),
	.z_real(y1_r),
	.z_imag(y1_i)
	);

endmodule
