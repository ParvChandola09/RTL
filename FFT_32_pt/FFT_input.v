module FFT_input #
	(
	parameter N = 16,
	parameter Q = 8
	)
	(
	input  signed  [N-1:0]  x0,
	input  signed  [N-1:0]  x1,
	input  signed  [N-1:0]  x2,
	input  signed  [N-1:0]  x3,
	input  signed  [N-1:0]  x4,
	input  signed  [N-1:0]  x5,
	input  signed  [N-1:0]  x6,
	input  signed  [N-1:0]  x7,
	input  signed  [N-1:0]  x8,
	input  signed  [N-1:0]  x9,
	input  signed  [N-1:0]  x10,
	input  signed  [N-1:0]  x11,
	input  signed  [N-1:0]  x12,
	input  signed  [N-1:0]  x13,
	input  signed  [N-1:0]  x14,
	input  signed  [N-1:0]  x15,
	input  signed  [N-1:0]  x16,
	input  signed  [N-1:0]  x17,
	input  signed  [N-1:0]  x18,
	input  signed  [N-1:0]  x19,
	input  signed  [N-1:0]  x20,
	input  signed  [N-1:0]  x21,
	input  signed  [N-1:0]  x22,
	input  signed  [N-1:0]  x23,
	input  signed  [N-1:0]  x24,
	input  signed  [N-1:0]  x25,
	input  signed  [N-1:0]  x26,
	input  signed  [N-1:0]  x27,
	input  signed  [N-1:0]  x28,
	input  signed  [N-1:0]  x29,
	input  signed  [N-1:0]  x30,
	input  signed  [N-1:0]  x31,

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

//intermediate stages
wire signed [N-1:0]ST1_r[2*N-1:0];
wire signed [N-1:0]ST1_i[2*N-1:0];	

//16 2-point stage

//1
butterfly_2 #(.N(N),.Q(Q)) stage1
	(
	.x0_r(x0),
	.x0_i('0),
	.x1_r(x16),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[0]),
	.y0_i(ST1_i[0]),
	.y1_r(ST1_r[1]),
	.y1_i(ST1_i[1])
	);

//2
butterfly_2 #(.N(N),.Q(Q)) stage2
	(
	.x0_r(x8),
	.x0_i('0),
	.x1_r(x24),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[2]),
	.y0_i(ST1_i[2]),
	.y1_r(ST1_r[3]),
	.y1_i(ST1_i[3])
	);

//3
butterfly_2 #(.N(N),.Q(Q)) stage3
	(
	.x0_r(x4),
	.x0_i('0),
	.x1_r(x20),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[4]),
	.y0_i(ST1_i[4]),
	.y1_r(ST1_r[5]),
	.y1_i(ST1_i[5])
	);

//4
butterfly_2 #(.N(N),.Q(Q)) stage4
	(
	.x0_r(x12),
	.x0_i('0),
	.x1_r(x28),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[6]),
	.y0_i(ST1_i[6]),
	.y1_r(ST1_r[7]),
	.y1_i(ST1_i[7])
	);

//5
butterfly_2 #(.N(N),.Q(Q)) stage5
	(
	.x0_r(x2),
	.x0_i('0),
	.x1_r(x18),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[8]),
	.y0_i(ST1_i[8]),
	.y1_r(ST1_r[9]),
	.y1_i(ST1_i[9])
	);

//6
butterfly_2 #(.N(N),.Q(Q)) stage6
	(
	.x0_r(x10),
	.x0_i('0),
	.x1_r(x26),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[10]),
	.y0_i(ST1_i[10]),
	.y1_r(ST1_r[11]),
	.y1_i(ST1_i[11])
	);

//7
butterfly_2 #(.N(N),.Q(Q)) stage7
	(
	.x0_r(x6),
	.x0_i('0),
	.x1_r(x22),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[12]),
	.y0_i(ST1_i[12]),
	.y1_r(ST1_r[13]),
	.y1_i(ST1_i[13])
	);

//8
butterfly_2 #(.N(N),.Q(Q)) stage8
	(
	.x0_r(x14),
	.x0_i('0),
	.x1_r(x30),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[14]),
	.y0_i(ST1_i[14]),
	.y1_r(ST1_r[15]),
	.y1_i(ST1_i[15])
	);

//9
butterfly_2 #(.N(N),.Q(Q)) stage9
	(
	.x0_r(x1),
	.x0_i('0),
	.x1_r(x17),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[16]),
	.y0_i(ST1_i[16]),
	.y1_r(ST1_r[17]),
	.y1_i(ST1_i[17])
	);

//10
butterfly_2 #(.N(N),.Q(Q)) stage10
	(
	.x0_r(x9),
	.x0_i('0),
	.x1_r(x25),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[18]),
	.y0_i(ST1_i[18]),
	.y1_r(ST1_r[19]),
	.y1_i(ST1_i[19])
	);

//11
butterfly_2 #(.N(N),.Q(Q)) stage11
	(
	.x0_r(x5),
	.x0_i('0),
	.x1_r(x21),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[20]),
	.y0_i(ST1_i[20]),
	.y1_r(ST1_r[21]),
	.y1_i(ST1_i[21])
	);

//12
butterfly_2 #(.N(N),.Q(Q)) stage12
	(
	.x0_r(x13),
	.x0_i('0),
	.x1_r(x29),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[22]),
	.y0_i(ST1_i[22]),
	.y1_r(ST1_r[23]),
	.y1_i(ST1_i[23])
	);

//13
butterfly_2 #(.N(N),.Q(Q)) stage13
	(
	.x0_r(x3),
	.x0_i('0),
	.x1_r(x19),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[24]),
	.y0_i(ST1_i[24]),
	.y1_r(ST1_r[25]),
	.y1_i(ST1_i[25])
	);

//14
butterfly_2 #(.N(N),.Q(Q)) stage14
	(
	.x0_r(x11),
	.x0_i('0),
	.x1_r(x27),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[26]),
	.y0_i(ST1_i[26]),
	.y1_r(ST1_r[27]),
	.y1_i(ST1_i[27])
	);

//15
butterfly_2 #(.N(N),.Q(Q)) stage15
	(
	.x0_r(x7),
	.x0_i('0),
	.x1_r(x23),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[28]),
	.y0_i(ST1_i[28]),
	.y1_r(ST1_r[29]),
	.y1_i(ST1_i[29])
	);

//16
butterfly_2 #(.N(N),.Q(Q)) stage16
	(
	.x0_r(x15),
	.x0_i('0),
	.x1_r(x31),
	.x1_i('0),
	.tw_r(16'b0000000100000000),
	.tw_i(16'b0000000000000000),
	.y0_r(ST1_r[30]),
	.y0_i(ST1_i[30]),
	.y1_r(ST1_r[31]),
	.y1_i(ST1_i[31])
	);

//2nd stage 8 4 point DFT
//1
wire signed [N-1:0]ST2_r[2*N-1:0];
wire signed [N-1:0]ST2_i[2*N-1:0];

butterfly_4 #(.N(N),.Q(Q)) st2_1
	(
	.x0_r(ST1_r[0]),
	.x0_i(ST1_i[0]),
	.x1_r(ST1_r[1]),
	.x1_i(ST1_i[1]),
	.x2_r(ST1_r[2]),
	.x2_i(ST1_i[2]),
	.x3_r(ST1_r[3]),
	.x3_i(ST1_i[3]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[0]),
	.y0_i(ST2_i[0]),
	.y1_r(ST2_r[1]),
	.y1_i(ST2_i[1]),
	.y2_r(ST2_r[2]),
	.y2_i(ST2_i[2]),
	.y3_r(ST2_r[3]),
	.y3_i(ST2_i[3])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_2
	(
	.x0_r(ST1_r[4]),
	.x0_i(ST1_i[4]),
	.x1_r(ST1_r[5]),
	.x1_i(ST1_i[5]),
	.x2_r(ST1_r[6]),
	.x2_i(ST1_i[6]),
	.x3_r(ST1_r[7]),
	.x3_i(ST1_i[7]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[4]),
	.y0_i(ST2_i[4]),
	.y1_r(ST2_r[5]),
	.y1_i(ST2_i[5]),
	.y2_r(ST2_r[6]),
	.y2_i(ST2_i[6]),
	.y3_r(ST2_r[7]),
	.y3_i(ST2_i[7])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_3
	(
	.x0_r(ST1_r[8]),
	.x0_i(ST1_i[8]),
	.x1_r(ST1_r[9]),
	.x1_i(ST1_i[9]),
	.x2_r(ST1_r[10]),
	.x2_i(ST1_i[10]),
	.x3_r(ST1_r[11]),
	.x3_i(ST1_i[11]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[8]),
	.y0_i(ST2_i[8]),
	.y1_r(ST2_r[9]),
	.y1_i(ST2_i[9]),
	.y2_r(ST2_r[10]),
	.y2_i(ST2_i[10]),
	.y3_r(ST2_r[11]),
	.y3_i(ST2_i[11])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_4
	(
	.x0_r(ST1_r[12]),
	.x0_i(ST1_i[12]),
	.x1_r(ST1_r[13]),
	.x1_i(ST1_i[13]),
	.x2_r(ST1_r[14]),
	.x2_i(ST1_i[14]),
	.x3_r(ST1_r[15]),
	.x3_i(ST1_i[15]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[12]),
	.y0_i(ST2_i[12]),
	.y1_r(ST2_r[13]),
	.y1_i(ST2_i[13]),
	.y2_r(ST2_r[14]),
	.y2_i(ST2_i[14]),
	.y3_r(ST2_r[15]),
	.y3_i(ST2_i[15])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_5
	(
	.x0_r(ST1_r[16]),
	.x0_i(ST1_i[16]),
	.x1_r(ST1_r[17]),
	.x1_i(ST1_i[17]),
	.x2_r(ST1_r[18]),
	.x2_i(ST1_i[18]),
	.x3_r(ST1_r[19]),
	.x3_i(ST1_i[19]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[16]),
	.y0_i(ST2_i[16]),
	.y1_r(ST2_r[17]),
	.y1_i(ST2_i[17]),
	.y2_r(ST2_r[18]),
	.y2_i(ST2_i[18]),
	.y3_r(ST2_r[19]),
	.y3_i(ST2_i[19])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_6
	(
	.x0_r(ST1_r[20]),
	.x0_i(ST1_i[20]),
	.x1_r(ST1_r[21]),
	.x1_i(ST1_i[21]),
	.x2_r(ST1_r[22]),
	.x2_i(ST1_i[22]),
	.x3_r(ST1_r[23]),
	.x3_i(ST1_i[23]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[20]),
	.y0_i(ST2_i[20]),
	.y1_r(ST2_r[21]),
	.y1_i(ST2_i[21]),
	.y2_r(ST2_r[22]),
	.y2_i(ST2_i[22]),
	.y3_r(ST2_r[23]),
	.y3_i(ST2_i[23])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_7
	(
	.x0_r(ST1_r[24]),
	.x0_i(ST1_i[24]),
	.x1_r(ST1_r[25]),
	.x1_i(ST1_i[25]),
	.x2_r(ST1_r[26]),
	.x2_i(ST1_i[26]),
	.x3_r(ST1_r[27]),
	.x3_i(ST1_i[27]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[24]),
	.y0_i(ST2_i[24]),
	.y1_r(ST2_r[25]),
	.y1_i(ST2_i[25]),
	.y2_r(ST2_r[26]),
	.y2_i(ST2_i[26]),
	.y3_r(ST2_r[27]),
	.y3_i(ST2_i[27])
	);

butterfly_4 #(.N(N),.Q(Q)) st2_8
	(
	.x0_r(ST1_r[28]),
	.x0_i(ST1_i[28]),
	.x1_r(ST1_r[29]),
	.x1_i(ST1_i[29]),
	.x2_r(ST1_r[30]),
	.x2_i(ST1_i[30]),
	.x3_r(ST1_r[31]),
	.x3_i(ST1_i[31]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000000000000),
	.tw2_i(16'b1111111100000000),
	.y0_r(ST2_r[28]),
	.y0_i(ST2_i[28]),
	.y1_r(ST2_r[29]),
	.y1_i(ST2_i[29]),
	.y2_r(ST2_r[30]),
	.y2_i(ST2_i[30]),
	.y3_r(ST2_r[31]),
	.y3_i(ST2_i[31])
	);


//3rd stage
wire signed [N-1:0]ST3_r[2*N-1:0];
wire signed [N-1:0]ST3_i[2*N-1:0];

butterfly_8 #(.N(N),.Q(Q)) st3_1
	(
	.x0_r(ST2_r[0]),
	.x0_i(ST2_i[0]),
	.x1_r(ST2_r[1]),
	.x1_i(ST2_i[1]),
	.x2_r(ST2_r[2]),
	.x2_i(ST2_i[2]),
	.x3_r(ST2_r[3]),
	.x3_i(ST2_i[3]),
	.x4_r(ST2_r[4]),
	.x4_i(ST2_i[4]),
	.x5_r(ST2_r[5]),
	.x5_i(ST2_i[5]),
	.x6_r(ST2_r[6]),
	.x6_i(ST2_i[6]),
	.x7_r(ST2_r[7]),
	.x7_i(ST2_i[7]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000010110101),
	.tw2_i(16'b1111111101001011),
	.tw3_r(16'b0000000000000000),
	.tw3_i(16'b1111111100000000),
	.tw4_r(16'b1111111101001011),
	.tw4_i(16'b1111111101001011),
	.y0_r(ST3_r[0]),
	.y0_i(ST3_i[0]),
	.y1_r(ST3_r[1]),
	.y1_i(ST3_i[1]),
	.y2_r(ST3_r[2]),
	.y2_i(ST3_i[2]),
	.y3_r(ST3_r[3]),
	.y3_i(ST3_i[3]),
	.y4_r(ST3_r[4]),
	.y4_i(ST3_i[4]),
	.y5_r(ST3_r[5]),
	.y5_i(ST3_i[5]),
	.y6_r(ST3_r[6]),
	.y6_i(ST3_i[6]),
	.y7_r(ST3_r[7]),
	.y7_i(ST3_i[7])
	);

butterfly_8 #(.N(N),.Q(Q)) st3_2
	(
	.x0_r(ST2_r[8]),
	.x0_i(ST2_i[8]),
	.x1_r(ST2_r[9]),
	.x1_i(ST2_i[9]),
	.x2_r(ST2_r[10]),
	.x2_i(ST2_i[10]),
	.x3_r(ST2_r[11]),
	.x3_i(ST2_i[11]),
	.x4_r(ST2_r[12]),
	.x4_i(ST2_i[12]),
	.x5_r(ST2_r[13]),
	.x5_i(ST2_i[13]),
	.x6_r(ST2_r[14]),
	.x6_i(ST2_i[14]),
	.x7_r(ST2_r[15]),
	.x7_i(ST2_i[15]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000010110101),
	.tw2_i(16'b1111111101001011),
	.tw3_r(16'b0000000000000000),
	.tw3_i(16'b1111111100000000),
	.tw4_r(16'b1111111101001011),
	.tw4_i(16'b1111111101001011),
	.y0_r(ST3_r[8]),
	.y0_i(ST3_i[8]),
	.y1_r(ST3_r[9]),
	.y1_i(ST3_i[9]),
	.y2_r(ST3_r[10]),
	.y2_i(ST3_i[10]),
	.y3_r(ST3_r[11]),
	.y3_i(ST3_i[11]),
	.y4_r(ST3_r[12]),
	.y4_i(ST3_i[12]),
	.y5_r(ST3_r[13]),
	.y5_i(ST3_i[13]),
	.y6_r(ST3_r[14]),
	.y6_i(ST3_i[14]),
	.y7_r(ST3_r[15]),
	.y7_i(ST3_i[15])
	);

butterfly_8 #(.N(N),.Q(Q)) st3_3
	(
	.x0_r(ST2_r[16]),
	.x0_i(ST2_i[16]),
	.x1_r(ST2_r[17]),
	.x1_i(ST2_i[17]),
	.x2_r(ST2_r[18]),
	.x2_i(ST2_i[18]),
	.x3_r(ST2_r[19]),
	.x3_i(ST2_i[19]),
	.x4_r(ST2_r[20]),
	.x4_i(ST2_i[20]),
	.x5_r(ST2_r[21]),
	.x5_i(ST2_i[21]),
	.x6_r(ST2_r[22]),
	.x6_i(ST2_i[22]),
	.x7_r(ST2_r[23]),
	.x7_i(ST2_i[23]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000010110101),
	.tw2_i(16'b1111111101001011),
	.tw3_r(16'b0000000000000000),
	.tw3_i(16'b1111111100000000),
	.tw4_r(16'b1111111101001011),
	.tw4_i(16'b1111111101001011),
	.y0_r(ST3_r[16]),
	.y0_i(ST3_i[16]),
	.y1_r(ST3_r[17]),
	.y1_i(ST3_i[17]),
	.y2_r(ST3_r[18]),
	.y2_i(ST3_i[18]),
	.y3_r(ST3_r[19]),
	.y3_i(ST3_i[19]),
	.y4_r(ST3_r[20]),
	.y4_i(ST3_i[20]),
	.y5_r(ST3_r[21]),
	.y5_i(ST3_i[21]),
	.y6_r(ST3_r[22]),
	.y6_i(ST3_i[22]),
	.y7_r(ST3_r[23]),
	.y7_i(ST3_i[23])
	);

butterfly_8 #(.N(N),.Q(Q)) st3_4
	(
	.x0_r(ST2_r[24]),
	.x0_i(ST2_i[24]),
	.x1_r(ST2_r[25]),
	.x1_i(ST2_i[25]),
	.x2_r(ST2_r[26]),
	.x2_i(ST2_i[26]),
	.x3_r(ST2_r[27]),
	.x3_i(ST2_i[27]),
	.x4_r(ST2_r[28]),
	.x4_i(ST2_i[28]),
	.x5_r(ST2_r[29]),
	.x5_i(ST2_i[29]),
	.x6_r(ST2_r[30]),
	.x6_i(ST2_i[30]),
	.x7_r(ST2_r[31]),
	.x7_i(ST2_i[31]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000010110101),
	.tw2_i(16'b1111111101001011),
	.tw3_r(16'b0000000000000000),
	.tw3_i(16'b1111111100000000),
	.tw4_r(16'b1111111101001011),
	.tw4_i(16'b1111111101001011),
	.y0_r(ST3_r[24]),
	.y0_i(ST3_i[24]),
	.y1_r(ST3_r[25]),
	.y1_i(ST3_i[25]),
	.y2_r(ST3_r[26]),
	.y2_i(ST3_i[26]),
	.y3_r(ST3_r[27]),
	.y3_i(ST3_i[27]),
	.y4_r(ST3_r[28]),
	.y4_i(ST3_i[28]),
	.y5_r(ST3_r[29]),
	.y5_i(ST3_i[29]),
	.y6_r(ST3_r[30]),
	.y6_i(ST3_i[30]),
	.y7_r(ST3_r[31]),
	.y7_i(ST3_i[31])
	);

//4th stage
wire signed [N-1:0]ST4_r[2*N-1:0];
wire signed [N-1:0]ST4_i[2*N-1:0];

butterfly_16 #(.N(N),.Q(Q)) st4_1
	(
	.x0_r(ST3_r[0]),
	.x0_i(ST3_i[0]),
	.x1_r(ST3_r[1]),
	.x1_i(ST3_i[1]),
	.x2_r(ST3_r[2]),
	.x2_i(ST3_i[2]),
	.x3_r(ST3_r[3]),
	.x3_i(ST3_i[3]),
	.x4_r(ST3_r[4]),
	.x4_i(ST3_i[4]),
	.x5_r(ST3_r[5]),
	.x5_i(ST3_i[5]),
	.x6_r(ST3_r[6]),
	.x6_i(ST3_i[6]),
	.x7_r(ST3_r[7]),
	.x7_i(ST3_i[7]),
	.x8_r(ST3_r[8]),
	.x8_i(ST3_i[8]),
	.x9_r(ST3_r[9]),
	.x9_i(ST3_i[9]),
	.x10_r(ST3_r[10]),
	.x10_i(ST3_i[10]),
	.x11_r(ST3_r[11]),
	.x11_i(ST3_i[11]),
	.x12_r(ST3_r[12]),
	.x12_i(ST3_i[12]),
	.x13_r(ST3_r[13]),
	.x13_i(ST3_i[13]),
	.x14_r(ST3_r[14]),
	.x14_i(ST3_i[14]),
	.x15_r(ST3_r[15]),
	.x15_i(ST3_i[15]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000011101101),
	.tw2_i(16'b1111111110011110),
	.tw3_r(16'b0000000010110101),
	.tw3_i(16'b1111111101001011),
	.tw4_r(16'b0000000001100010),
	.tw4_i(16'b1111111100010011),
	.tw5_r(16'b0000000000000000),
	.tw5_i(16'b1111111100000000),
	.tw6_r(16'b1111111110011110),
	.tw6_i(16'b1111111100010011),
	.tw7_r(16'b1111111101001011),
	.tw7_i(16'b1111111101001011),
	.tw8_r(16'b1111111100010011),
	.tw8_i(16'b1111111110011110),
	.y0_r(ST4_r[0]),
	.y0_i(ST4_i[0]),
	.y1_r(ST4_r[1]),
	.y1_i(ST4_i[1]),
	.y2_r(ST4_r[2]),
	.y2_i(ST4_i[2]),
	.y3_r(ST4_r[3]),
	.y3_i(ST4_i[3]),
	.y4_r(ST4_r[4]),
	.y4_i(ST4_i[4]),
	.y5_r(ST4_r[5]),
	.y5_i(ST4_i[5]),
	.y6_r(ST4_r[6]),
	.y6_i(ST4_i[6]),
	.y7_r(ST4_r[7]),
	.y7_i(ST4_i[7]),
	.y8_r(ST4_r[8]),
	.y8_i(ST4_i[8]),
	.y9_r(ST4_r[9]),
	.y9_i(ST4_i[9]),
	.y10_r(ST4_r[10]),
	.y10_i(ST4_i[10]),
	.y11_r(ST4_r[11]),
	.y11_i(ST4_i[11]),
	.y12_r(ST4_r[12]),
	.y12_i(ST4_i[12]),
	.y13_r(ST4_r[13]),
	.y13_i(ST4_i[13]),
	.y14_r(ST4_r[14]),
	.y14_i(ST4_i[14]),
	.y15_r(ST4_r[15]),
	.y15_i(ST4_i[15])
	);

butterfly_16 #(.N(N),.Q(Q)) st4_2
	(
	.x0_r(ST3_r[16]),
	.x0_i(ST3_i[16]),
	.x1_r(ST3_r[17]),
	.x1_i(ST3_i[17]),
	.x2_r(ST3_r[18]),
	.x2_i(ST3_i[18]),
	.x3_r(ST3_r[19]),
	.x3_i(ST3_i[19]),
	.x4_r(ST3_r[20]),
	.x4_i(ST3_i[20]),
	.x5_r(ST3_r[21]),
	.x5_i(ST3_i[21]),
	.x6_r(ST3_r[22]),
	.x6_i(ST3_i[22]),
	.x7_r(ST3_r[23]),
	.x7_i(ST3_i[23]),
	.x8_r(ST3_r[24]),
	.x8_i(ST3_i[24]),
	.x9_r(ST3_r[25]),
	.x9_i(ST3_i[25]),
	.x10_r(ST3_r[26]),
	.x10_i(ST3_i[26]),
	.x11_r(ST3_r[27]),
	.x11_i(ST3_i[27]),
	.x12_r(ST3_r[28]),
	.x12_i(ST3_i[28]),
	.x13_r(ST3_r[29]),
	.x13_i(ST3_i[29]),
	.x14_r(ST3_r[30]),
	.x14_i(ST3_i[30]),
	.x15_r(ST3_r[31]),
	.x15_i(ST3_i[31]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000011101101),
	.tw2_i(16'b1111111110011110),
	.tw3_r(16'b0000000010110101),
	.tw3_i(16'b1111111101001011),
	.tw4_r(16'b0000000001100010),
	.tw4_i(16'b1111111100010011),
	.tw5_r(16'b0000000000000000),
	.tw5_i(16'b1111111100000000),
	.tw6_r(16'b1111111110011110),
	.tw6_i(16'b1111111100010011),
	.tw7_r(16'b1111111101001011),
	.tw7_i(16'b1111111101001011),
	.tw8_r(16'b1111111100010011),
	.tw8_i(16'b1111111110011110),
	.y0_r(ST4_r[16]),
	.y0_i(ST4_i[16]),
	.y1_r(ST4_r[17]),
	.y1_i(ST4_i[17]),
	.y2_r(ST4_r[18]),
	.y2_i(ST4_i[18]),
	.y3_r(ST4_r[19]),
	.y3_i(ST4_i[19]),
	.y4_r(ST4_r[20]),
	.y4_i(ST4_i[20]),
	.y5_r(ST4_r[21]),
	.y5_i(ST4_i[21]),
	.y6_r(ST4_r[22]),
	.y6_i(ST4_i[22]),
	.y7_r(ST4_r[23]),
	.y7_i(ST4_i[23]),
	.y8_r(ST4_r[24]),
	.y8_i(ST4_i[24]),
	.y9_r(ST4_r[25]),
	.y9_i(ST4_i[25]),
	.y10_r(ST4_r[26]),
	.y10_i(ST4_i[26]),
	.y11_r(ST4_r[27]),
	.y11_i(ST4_i[27]),
	.y12_r(ST4_r[28]),
	.y12_i(ST4_i[28]),
	.y13_r(ST4_r[29]),
	.y13_i(ST4_i[29]),
	.y14_r(ST4_r[30]),
	.y14_i(ST4_i[30]),
	.y15_r(ST4_r[31]),
	.y15_i(ST4_i[31])
	);


//5th stage
butterfly_32 #(.N(N),.Q(Q)) st5_1
	(
	.x0_r(ST4_r[0]),
	.x0_i(ST4_i[0]),
	.x1_r(ST4_r[1]),
	.x1_i(ST4_i[1]),
	.x2_r(ST4_r[2]),
	.x2_i(ST4_i[2]),
	.x3_r(ST4_r[3]),
	.x3_i(ST4_i[3]),
	.x4_r(ST4_r[4]),
	.x4_i(ST4_i[4]),
	.x5_r(ST4_r[5]),
	.x5_i(ST4_i[5]),
	.x6_r(ST4_r[6]),
	.x6_i(ST4_i[6]),
	.x7_r(ST4_r[7]),
	.x7_i(ST4_i[7]),
	.x8_r(ST4_r[8]),
	.x8_i(ST4_i[8]),
	.x9_r(ST4_r[9]),
	.x9_i(ST4_i[9]),
	.x10_r(ST4_r[10]),
	.x10_i(ST4_i[10]),
	.x11_r(ST4_r[11]),
	.x11_i(ST4_i[11]),
	.x12_r(ST4_r[12]),
	.x12_i(ST4_i[12]),
	.x13_r(ST4_r[13]),
	.x13_i(ST4_i[13]),
	.x14_r(ST4_r[14]),
	.x14_i(ST4_i[14]),
	.x15_r(ST4_r[15]),
	.x15_i(ST4_i[15]),
	.x16_r(ST4_r[16]),
	.x16_i(ST4_i[16]),
	.x17_r(ST4_r[17]),
	.x17_i(ST4_i[17]),
	.x18_r(ST4_r[18]),
	.x18_i(ST4_i[18]),
	.x19_r(ST4_r[19]),
	.x19_i(ST4_i[19]),
	.x20_r(ST4_r[20]),
	.x20_i(ST4_i[20]),
	.x21_r(ST4_r[21]),
	.x21_i(ST4_i[21]),
	.x22_r(ST4_r[22]),
	.x22_i(ST4_i[22]),
	.x23_r(ST4_r[23]),
	.x23_i(ST4_i[23]),
	.x24_r(ST4_r[24]),
	.x24_i(ST4_i[24]),
	.x25_r(ST4_r[25]),
	.x25_i(ST4_i[25]),
	.x26_r(ST4_r[26]),
	.x26_i(ST4_i[26]),
	.x27_r(ST4_r[27]),
	.x27_i(ST4_i[27]),
	.x28_r(ST4_r[28]),
	.x28_i(ST4_i[28]),
	.x29_r(ST4_r[29]),
	.x29_i(ST4_i[29]),
	.x30_r(ST4_r[30]),
	.x30_i(ST4_i[30]),
	.x31_r(ST4_r[31]),
	.x31_i(ST4_i[31]),
	.tw1_r(16'b0000000100000000),
	.tw1_i(16'b0000000000000000),
	.tw2_r(16'b0000000011111011),
	.tw2_i(16'b1111111111001110),
	.tw3_r(16'b0000000011101101),
	.tw3_i(16'b1111111110011110),
	.tw4_r(16'b0000000011010101),
	.tw4_i(16'b1111111101110010),
	.tw5_r(16'b0000000010110101),
	.tw5_i(16'b1111111101001011),
	.tw6_r(16'b0000000010001110),
	.tw6_i(16'b1111111100101011),
	.tw7_r(16'b0000000001100010),
	.tw7_i(16'b1111111100010011),
	.tw8_r(16'b0000000000110010),
	.tw8_i(16'b1111111100000101),
	.tw9_r(16'b0000000000000000),
	.tw9_i(16'b1111111100000000),
	.tw10_r(16'b1111111111001110),
	.tw10_i(16'b1111111100000101),
	.tw11_r(16'b1111111110011110),
	.tw11_i(16'b1111111100010011),
	.tw12_r(16'b1111111101110010),
	.tw12_i(16'b1111111100101011),
	.tw13_r(16'b1111111101001011),
	.tw13_i(16'b1111111101001011),
	.tw14_r(16'b1111111100101011),
	.tw14_i(16'b1111111101110010),
	.tw15_r(16'b1111111100010011),
	.tw15_i(16'b1111111110011110),
	.tw16_r(16'b1111111100000101),
	.tw16_i(16'b1111111111001110),
	.y0_r(y0_r),
	.y0_i(y0_i),
	.y1_r(y1_r),
	.y1_i(y1_i),
	.y2_r(y2_r),
	.y2_i(y2_i),
	.y3_r(y3_r),
	.y3_i(y3_i),
	.y4_r(y4_r),
	.y4_i(y4_i),
	.y5_r(y5_r),
	.y5_i(y5_i),
	.y6_r(y6_r),
	.y6_i(y6_i),
	.y7_r(y7_r),
	.y7_i(y7_i),
	.y8_r(y8_r),
	.y8_i(y8_i),
	.y9_r(y9_r),
	.y9_i(y9_i),
	.y10_r(y10_r),
	.y10_i(y10_i),
	.y11_r(y11_r),
	.y11_i(y11_i),
	.y12_r(y12_r),
	.y12_i(y12_i),
	.y13_r(y13_r),
	.y13_i(y13_i),
	.y14_r(y14_r),
	.y14_i(y14_i),
	.y15_r(y15_r),
	.y15_i(y15_i),
	.y16_r(y16_r),
	.y16_i(y16_i),
	.y17_r(y17_r),
	.y17_i(y17_i),
	.y18_r(y18_r),
	.y18_i(y18_i),
	.y19_r(y19_r),
	.y19_i(y19_i),
	.y20_r(y20_r),
	.y20_i(y20_i),
	.y21_r(y21_r),
	.y21_i(y21_i),
	.y22_r(y22_r),
	.y22_i(y22_i),
	.y23_r(y23_r),
	.y23_i(y23_i),
	.y24_r(y24_r),
	.y24_i(y24_i),
	.y25_r(y25_r),
	.y25_i(y25_i),
	.y26_r(y26_r),
	.y26_i(y26_i),
	.y27_r(y27_r),
	.y27_i(y27_i),
	.y28_r(y28_r),
	.y28_i(y28_i),
	.y29_r(y29_r),
	.y29_i(y29_i),
	.y30_r(y30_r),
	.y30_i(y30_i),
	.y31_r(y31_r),
	.y31_i(y31_i)
	);

endmodule



////////////////////////////////////
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


///////////////////////////////////////
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

///////////////////////////////////////////
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

////////////////////////////

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

////////////////////////////
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
