module complex_subtractor#
	(
	parameter N = 16
	)
	(
	input signed  [N-1:0] x_real,
	input signed  [N-1:0] x_imag,
	input signed  [N-1:0] y_real,
	input signed  [N-1:0] y_imag,
	output reg signed [N-1:0] z_real,
	output reg signed [N-1:0] z_imag
	);


assign z_real = x_real - y_real;
assign z_imag = x_imag - y_imag;

endmodule
