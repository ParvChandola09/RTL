/*#####################################################################################
#######################################################################################
		    VERILOG PROGRAM FOR FIR FILTER DESIGN
#######################################################################################
#######################################################################################
1.Filter coefficients and the inputs are generated in MATLAB.
2.Filter specifications passband=100Hz, Stopband = 400Hz, Order = 15
3.Input applied is a combination of two sinusoidal signals
  x1 = 200Hz and x2 = 600Hz
4.Output of the filter only contains the signal x1 while x2 is filtered out.
5.In the current code 15 bits are fractional bits and 16th is sign bit.(Fixed point annotation)
 */

module test_FIR(clk, reset, data_in, data_out);

parameter N = 16;
parameter M = 32;

input clk, reset;
input signed [N-1:0] data_in;
output reg signed [N-1:0] data_out;

//filter coefficients generated through matlab and converted to signed binary.

wire signed [N-1:0]h0 = 16'b1111111110001111;
wire signed [N-1:0]h1 = 16'b1111110111101001;
wire signed [N-1:0]h2 = 16'b1111110011111001;
wire signed [N-1:0]h3 = 16'b1111111000110000;
wire signed [N-1:0]h4 = 16'b0000001110101110;
wire signed [N-1:0]h5 = 16'b0000110100101001;
wire signed [N-1:0]h6 = 16'b0001011110100101;
wire signed [N-1:0]h7 = 16'b0001111010010100;
wire signed [N-1:0]h8 = 16'b0001111010010100;
wire signed [N-1:0]h9 = 16'b0001011110100101;
wire signed [N-1:0]h10 = 16'b0000110100101001;
wire signed [N-1:0]h11 = 16'b0000001110101110;
wire signed [N-1:0]h12 = 16'b1111111000110000;
wire signed [N-1:0]h13 = 16'b1111110011111001;
wire signed [N-1:0]h14 = 16'b1111110111101001;
wire signed [N-1:0]h15 = 16'b1111111110001111;


reg signed [M-1:0] M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15; //for storing the multiplication results
reg signed [N-1:0] x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15; // buffer resistors for x[n-1],x[n-2]...etc
reg signed [M-1:0] y; //storing the output after accumulation

always @(posedge clk)
begin

if(reset)
begin
M0  <= 32'd0;
M1  <= 32'd0;
M2  <= 32'd0;
M3  <= 32'd0;
M4  <= 32'd0;
M5  <= 32'd0;
M6  <= 32'd0;
M7  <= 32'd0;
M8  <= 32'd0;
M9  <= 32'd0;
M10 <= 32'd0;
M11 <= 32'd0;
M12 <= 32'd0;
M13 <= 32'd0;
M14 <= 32'd0;
M15 <= 32'd0;

y <= 32'd0;

x0  <= 16'd0;
x1  <= 16'd0;
x2  <= 16'd0;
x3  <= 16'd0;
x4  <= 16'd0;
x5  <= 16'd0;
x6  <= 16'd0;
x7  <= 16'd0;
x8  <= 16'd0;
x9  <= 16'd0;
x10 <= 16'd0;
x11 <= 16'd0;
x12 <= 16'd0;
x13 <= 16'd0;
x14 <= 16'd0;
x15 <= 16'd0;
end

else
begin
//shifiting the values delays
x0  <= data_in;
x1  <= x0;  //x[n-1]
x2  <= x1;  //x[n-2]
x3  <= x2;  //x[n-3]
x4  <= x3;
x5  <= x4;
x6  <= x5;
x7  <= x6;
x8  <= x7;
x9  <= x8;
x10 <= x9;
x11 <= x10;
x12 <= x11;
x13 <= x12;
x14 <= x13;
x15 <= x14;

//The multiplication block
M0   <= x0*h0;
M1   <= x1*h1;
M2   <= x2*h2;
M3   <= x3*h3;
M4   <= x4*h4;
M5   <= x5*h5;
M6   <= x6*h6;
M7   <= x7*h7;
M8   <= x8*h8;
M9   <= x9*h9;
M10  <= x10*h10;
M11  <= x11*h11;
M12  <= x12*h12;
M13  <= x13*h13;
M14  <= x14*h14;
M15  <= x15*h15;

//accumulation block

y <= M0+M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11+M12+M13+M14+M15;

end

end
assign data_out = y[30:15]; //truncating the output bits.
endmodule


//TESTBENCH for the FIR filter

module test_FIR_tb;

parameter N=16;
parameter M=32;

reg signed [N-1:0]data_in;
wire signed [N-1:0]data_out;
reg clk,reset;
reg [11:0]Address;
reg signed [15:0]RAMM[2000:0];
integer mcd1;

test_FIR DUT(clk,reset,data_in,data_out);

//reading the input file
initial
begin
$readmemb("/home/annu/Parv/verilog_progs/FIR_filter/input_comb.dat", RAMM);
end

//writing on the file
initial begin
mcd1 = $fopen("filter_output_comb1.txt","w");
end

initial
begin
clk = 1'b0;
end

always #10 clk <= ~clk;

initial
begin
reset = 1'b1; 
#15 reset = 1'b0; data_in = 0;
#15 Address = 0;
end

always@ (posedge clk) begin
data_in <= RAMM[Address];
if(data_out==16'hxxxx)
$fwriteb(mcd1,"%d \n",0);
else
$fwriteb(mcd1,"%d \n",data_out);
end

//initial
//Address = 1;
always @(posedge clk)
begin
if(Address == 2000)
Address = 0;
else
Address = Address+1;
end

endmodule