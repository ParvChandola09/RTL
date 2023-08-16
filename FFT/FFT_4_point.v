//8-point DFT calcultion
module DFT8(
	input reg signed[15:0]x0,
	input reg signed[15:0]x1,
	input reg signed[15:0]x2,
	input reg signed[15:0]x3,
	input reg signed[15:0]x4,
	input reg signed[15:0]x5,
	input reg signed[15:0]x6,
	input reg signed[15:0]x7,
	output signed [15:0]X0_real,
	output signed [15:0]X0_imag,
	output signed [15:0]X1_real,
	output signed [15:0]X1_imag,
	output signed [15:0]X2_real,
	output signed [15:0]X2_imag,
	output signed [15:0]X3_real,
	output signed [15:0]X3_imag,
	output signed [15:0]X4_real,
	output signed [15:0]X4_imag,
	output signed [15:0]X5_real,
	output signed [15:0]X5_imag,
	output signed [15:0]X6_real,
	output signed [15:0]X6_imag,
	output signed [15:0]X7_real,
	output signed [15:0]X7_imag);

//diectly getting output of 2nd stage using 4 point FFT module
//upper box
wire signed[15:0]X20_real,X20_imag,X21_real,X21_imag,X22_real,X22_imag,X23_real,X23_imag;
DFT4 S21(x0,x2,x4,x6,X20_real,X20_imag,X21_real,X21_imag,X22_real,X22_imag,X23_real,X23_imag);

//lower box
wire signed[15:0]X24_real,X24_imag,X25_real,X25_imag,X26_real,X26_imag,X27_real,X27_imag;
DFT4 S22(x1,x3,x5,x7,X24_real,X24_imag,X25_real,X25_imag,X26_real,X26_imag,X27_real,X27_imag);

//twiddle factor caculations (8,8) fixed point representations
wire signed[15:0]W[5:0];
assign W[0] = 16'b0000000010110101;  //W_8^1 real part 
assign W[1] = 16'b1111111101001011;  //W_8^1 imag part
assign W[2] = 16'b0000000000000000;  //W_8^2 real part
assign W[3] = 16'b1111111100000000;  //W_8^2 imag part
assign W[4] = 16'b1111111101001011;  //W_8^3 real part
assign W[5] = 16'b1111111101001011;  //W_8^3 imag part

//3rd stage calculations
//X0
complex_adder CA31(X20_real,X20_imag,X24_real,X24_imag,X0_real,X0_imag);

//X1
wire signed[15:0]X25_r,X25_i;
complex_multiplier CA21(X25_real,X25_imag,W[0],W[1],X25_r,X25_i);
complex_adder CA32(X21_real,X21_imag,X25_r,X25_i,X1_real,X1_imag);

//X2
wire signed[15:0]X26_r,X26_i;
complex_multiplier CA22(X26_real,X26_imag,W[2],W[3],X26_r,X26_i);
complex_adder CA33(X22_real,X22_imag,X26_r,X26_i,X2_real,X2_imag);

//X3
wire signed[15:0]X27_r,X27_i;
complex_multiplier CA23(X27_real,X27_imag,W[4],W[5],X27_r,X27_i);
complex_adder CA34(X23_real,X23_imag,X27_r,X27_i,X3_real,X3_imag);

//X4
complex_subtractor CS31(X20_real,X20_imag,X24_real,X24_imag,X4_real,X4_imag);

//X5
complex_subtractor CS32(X21_real,X21_imag,X25_r,X25_i,X5_real,X5_imag);

//X6
complex_subtractor CS33(X22_real,X22_imag,X26_r,X26_i,X6_real,X6_imag);

//X7
complex_subtractor CS34(X23_real,X23_imag,X27_r,X27_i,X7_real,X7_imag);

endmodule



//4-point DFT calculation
module DFT4(
	input reg signed[15:0]x0,
	input reg signed[15:0]x1,
	input reg signed[15:0]x2,
	input reg signed[15:0]x3,
	output signed [15:0]X0_real,
	output signed [15:0]X0_imag,
	output signed [15:0]X1_real,
	output signed [15:0]X1_imag,
	output signed [15:0]X2_real,
	output signed [15:0]X2_imag,
	output signed [15:0]X3_real,
	output signed [15:0]X3_imag);

//Calculating the first stage output

//Upper box calculation
wire signed[15:0]X10_real,X10_imag,X11_real,X11_imag;
DFT2 ST11(x0,x2,X10_real,X10_imag,X11_real,X11_imag);

//Lower box calculation
wire signed[15:0]X12_real,X12_imag,X13_real,X13_imag;
DFT2 ST12(x1,x3,X12_real,X12_imag,X13_real,X13_imag);

//twiddle factor definition W_{8}^2
wire signed [15:0]W82_real = 16'b0000000000000000;
wire signed [15:0]W82_imag = 16'b1111111100000000;

//Calculation of 2nd stage
//X0
complex_adder CA21(X10_real,X10_imag,X12_real,X12_imag,X0_real,X0_imag);

//X1
wire signed[15:0]X13_r,X13_i;
complex_multiplier CM1(X13_real,X13_imag,W82_real,W82_imag,X13_r,X13_i);
complex_adder CA22(X11_real,X11_imag,X13_r,X13_i,X1_real,X1_imag);

//X2
complex_subtractor CS21(X10_real,X10_imag,X12_real,X12_imag,X2_real,X2_imag);

//X3
//using X13_r and X13_i already multiplication is done
complex_subtractor CS22(X11_real,X11_imag,X13_r,X13_i,X3_real,X3_imag);


endmodule


//2-point DFT 1st stage butterfly
module DFT2(
	input reg signed[15:0]x0,
	input reg signed[15:0]x1,
	output signed [15:0]X0_real,
	output signed [15:0]X0_imag,
	output signed [15:0]X1_real,
	output signed [15:0]X1_imag);

//X0
complex_adder CA11(x0,'0,x1,'0,X0_real,X0_imag);
//X1
complex_subtractor CS11(x0,'0,x1,'0,X1_real,X1_imag);


endmodule


//complex adder
module complex_adder(
	input reg signed[15:0]x_real,
	input reg signed[15:0]x_imag,
	input reg signed[15:0]y_real,
	input reg signed[15:0]y_imag,
	output reg signed[15:0]z_real,
	output reg signed[15:0]z_imag);

assign z_real = x_real+y_real;
assign z_imag = x_imag+y_imag;

endmodule


//complex subtractor
module complex_subtractor(
	input reg signed[15:0]x_real,
	input reg signed[15:0]x_imag,
	input reg signed[15:0]y_real,
	input reg signed[15:0]y_imag,
	output reg signed[15:0]z_real,
	output reg signed[15:0]z_imag);

assign z_real = x_real-y_real;
assign z_imag = x_imag-y_imag;

endmodule

//complex multipier of two f(8,8) fixed point multiplications
module complex_multiplier(
	input reg signed[15:0]x1_real,
	input reg signed[15:0]x1_imag,
	input reg signed[15:0]x2_real,
	input reg signed[15:0]x2_imag,
	output signed[15:0]z_real,
	output signed[15:0]z_imag);

wire signed[31:0]a,b,c,d,y1,y2;
assign a = x1_real*x2_real;
assign b = x1_imag*x2_imag;
assign c = x1_real*x2_imag;
assign d = x2_real*x1_imag;

assign y1 = a-b;
assign y2 = c+d;

assign z_real = y1[23:8];
assign z_imag = y2[23:8];

endmodule


module DFT8_tb();
reg signed[15:0]x[7:0];
wire signed[15:0]X[15:0];
integer fd;
integer i;
DFT8 DUT(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8],X[9],X[10],X[11],X[12],X[13],X[14],X[15]);
initial
begin
fd = $fopen("FFT8.txt","w");
x[0] = 16'b0000000100000000;
x[1] = 16'b1111111100000000;
x[2] = 16'b0000000100000000;
x[3] = 16'b1111111100000000;
x[4] = 16'b0000000000000000;
x[5] = 16'b0000000000000000;
x[6] = 16'b0000000000000000;
x[7] = 16'b0000000000000000;

#10
x[0] = 16'b0000000100000000;
x[1] = 16'b1111111100000000;
x[2] = 16'b0000000100000000;
x[3] = 16'b1111111100000000;
x[4] = 16'b0000000000000000;
x[5] = 16'b0000000000000000;
x[6] = 16'b0000000000000000;
x[7] = 16'b0000000000000000;
for(i=0;i<15;i=i+2)begin
	$fwriteb(fd,"%d   %d\n",X[i],X[i+1]);
end
#10 $fclose(fd);


end
endmodule

/*module DFT4_tb();
reg signed [15:0]x0;
reg signed [15:0]x1;
reg signed [15:0]x2;
reg signed [15:0]x3;
wire signed[15:0]X[7:0];
integer fd;
integer i;

DFT4 DUT(x0,x1,x2,x3,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7]);

initial
begin
fd = $fopen("FFT4.txt","w");
x0=4;x1=3;x2=2;x3=1;
#10 x0=4;x1=3;x2=2;x3=1;
for(i=0;i<7;i=i+2)begin
	$fwriteb(fd,"%d + j(%d)\n",X[i],X[i+1]);
end
#10 $fclose(fd);
end
endmodule*/

//testbench
/*module DFT2_tb();
reg signed [15:0]x0;
reg signed [15:0]x1;
wire signed[15:0]X[3:0];
integer fd;
integer i;
/*wire signed[16:0]X0_imag;
wire signed[16:0]X1_real;
wire signed[16:0]X1_imag;

DFT2 DUT(x0,x1,X[0],X[1],X[2],X[3]);

initial
begin
fd = $fopen("FFT2.txt", "w");
x0=2;x1=1;
#10 x0=2;x1=1;
for(i=0;i<3;i=i+2)begin
	$fwriteb(fd,"%d + j(%d)\n",X[i],X[i+1]);
end
#10 $fclose(fd);
end
endmodule*/






