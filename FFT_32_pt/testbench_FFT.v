module testbench_FFT();
parameter N = 16;
parameter Q = 8;
reg signed [N-1:0]x[2*N-1:0];
wire signed [N-1:0]Y_R[2*N-1:0];
wire signed [N-1:0]Y_I[2*N-1:0];
integer fd;
integer i;

FFT_input #(.N(N),.Q(Q)) hello
	(
	.x0(x[0]),
	.x1(x[1]),
	.x2(x[2]),
	.x3(x[3]),
	.x4(x[4]),
	.x5(x[5]),
	.x6(x[6]),
	.x7(x[7]),
	.x8(x[8]),
	.x9(x[9]),
	.x10(x[10]),
	.x11(x[11]),
	.x12(x[12]),
	.x13(x[13]),
	.x14(x[14]),
	.x15(x[15]),
	.x16(x[16]),
	.x17(x[17]),
	.x18(x[18]),
	.x19(x[19]),
	.x20(x[20]),
	.x21(x[21]),
	.x22(x[22]),
	.x23(x[23]),
	.x24(x[24]),
	.x25(x[25]),
	.x26(x[26]),
	.x27(x[27]),
	.x28(x[28]),
	.x29(x[29]),
	.x30(x[30]),
	.x31(x[31]),

	.y0_r(Y_R[0]),
	.y0_i(Y_I[0]),
	.y1_r(Y_R[1]),
	.y1_i(Y_I[1]),
	.y2_r(Y_R[2]),
	.y2_i(Y_I[2]),
	.y3_r(Y_R[3]),
	.y3_i(Y_I[3]),
	.y4_r(Y_R[4]),
	.y4_i(Y_I[4]),
	.y5_r(Y_R[5]),
	.y5_i(Y_I[5]),
	.y6_r(Y_R[6]),
	.y6_i(Y_I[6]),
	.y7_r(Y_R[7]),
	.y7_i(Y_I[7]),
	.y8_r(Y_R[8]),
	.y8_i(Y_I[8]),
	.y9_r(Y_R[9]),
	.y9_i(Y_I[9]),
	.y10_r(Y_R[10]),
	.y10_i(Y_I[10]),
	.y11_r(Y_R[11]),
	.y11_i(Y_I[11]),
	.y12_r(Y_R[12]),
	.y12_i(Y_I[12]),
	.y13_r(Y_R[13]),
	.y13_i(Y_I[13]),
	.y14_r(Y_R[14]),
	.y14_i(Y_I[14]),
	.y15_r(Y_R[15]),
	.y15_i(Y_I[15]),
	.y16_r(Y_R[16]),
	.y16_i(Y_I[16]),
	.y17_r(Y_R[17]),
	.y17_i(Y_I[17]),
	.y18_r(Y_R[18]),
	.y18_i(Y_I[18]),
	.y19_r(Y_R[19]),
	.y19_i(Y_I[19]),
	.y20_r(Y_R[20]),
	.y20_i(Y_I[20]),
	.y21_r(Y_R[21]),
	.y21_i(Y_I[21]),
	.y22_r(Y_R[22]),
	.y22_i(Y_I[22]),
	.y23_r(Y_R[23]),
	.y23_i(Y_I[23]),
	.y24_r(Y_R[24]),
	.y24_i(Y_I[24]),
	.y25_r(Y_R[25]),
	.y25_i(Y_I[25]),
	.y26_r(Y_R[26]),
	.y26_i(Y_I[26]),
	.y27_r(Y_R[27]),
	.y27_i(Y_I[27]),
	.y28_r(Y_R[28]),
	.y28_i(Y_I[28]),
	.y29_r(Y_R[29]),
	.y29_i(Y_I[29]),
	.y30_r(Y_R[30]),
	.y30_i(Y_I[30]),
	.y31_r(Y_R[31]),
	.y31_i(Y_I[31])
	);

initial
begin
$readmemb("/home/annu/Parv/verilog_progs/input_data.dat", x);
end

initial 
begin
fd = $fopen("FFT32.txt","w");
$readmemb("/home/annu/Parv/verilog_progs/input_data.dat", x);
#10
$readmemb("/home/annu/Parv/verilog_progs/input_data.dat", x);


for(i=0;i<32;i=i+1)begin
	$fwriteb(fd,"%d   %d\n",Y_R[i],Y_I[i]);
end
#10 $fclose(fd);

end


endmodule
