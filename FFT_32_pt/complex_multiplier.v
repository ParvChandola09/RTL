module complex_multiplier #
	(
	parameter N = 16,
	parameter Q = 8
	)
	(
	input signed  [N-1:0] x1_real,
	input signed  [N-1:0] x1_imag,
	input signed  [N-1:0] x2_real,
	input signed  [N-1:0] x2_imag,
	output reg signed [N-1:0] z_real,
	output reg signed [N-1:0] z_imag);

reg signed[2*N-1:0]a,b,c,d,y1,y2;

assign a = x1_real*x2_real;
assign b = x1_imag*x2_imag;
assign c = x1_real*x2_imag;
assign d = x2_real*x1_imag;

assign y1 = a-b;
assign y2 = c+d;

assign z_real = y1[23:8];
assign z_imag = y2[23:8];

endmodule
