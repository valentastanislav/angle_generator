float Fk(int rank, int L1, int L2, double spin1, double spin2){
// !   Returns values of standart function Fk which are present in the
// !   equation for angular correlation.
// !   Data are from some russian book.
// !   spin1 -- L1 --> spin2 -- L2 --> spin3
// !   Standard factors F2 are stored in the form F2[L1][L2][2*spin1][2*spin2]
// !   (reason is to include half-integer spins).
// !   At this moment only combinations (L=1,L'=1), (L=1,L'=2) and
// !   (L=2,L'=2) are assumed (no octupole mixing).
// !   F4 factors are stored only for L=2, L'=2 in the form F4[2*spin1][2*spin2].
  float  F2[5][5][25][21];
  float  F4[25][21];
  for(int k=0; k<25; ++k){
    for(int l=0; l<21; ++l){
      F4[k][l]=0.0;
      for(int i=0; i<5; ++i){
        for(int j=0; j<5; ++j){
          F2[i][j][k][l]=0.0;
        }
      }
    }
  }

  if (L1==0) {L1=1;}
  if (L2==0) {L2=1;}

  F2[1][1][0][2]=0.707;
  F2[1][1][2][2]=-0.354;
  F2[1][1][4][2]=0.071;
  F2[1][1][2][4]=0.418;
  F2[1][1][4][4]=-0.418;
  F2[1][1][6][4]=0.120;
  F2[1][1][4][6]=0.346;
  F2[1][1][6][6]=-0.433;
  F2[1][1][8][6]=0.144;
  F2[1][1][6][8]=0.313;
  F2[1][1][8][8]=-0.439;
  F2[1][1][10][8]=0.160;
  F2[1][1][8][10]=0.294;
  F2[1][1][10][10]=-0.442;
  F2[1][1][12][10]=0.170;
  F2[1][1][10][12]=0.282;
  F2[1][1][12][12]=-0.443;
  F2[1][1][14][12]=0.177;
  F2[1][1][12][14]=0.273;
  F2[1][1][14][14]=-0.444;
  F2[1][1][16][14]=0.183;
  F2[1][1][14][16]=0.267;
  F2[1][1][16][16]=-0.445;
  F2[1][1][18][16]=0.187;
  F2[1][1][16][18]=0.262;
  F2[1][1][18][18]=-0.445;
  F2[1][1][20][18]=0.191;
  F2[1][1][18][20]=0.258;
  F2[1][1][20][20]=-0.446;
  F2[1][1][22][20]=0.194;

  F2[1][2][2][2]=-1.061;
  F2[1][2][4][2]=0.474;
  F2[1][2][2][4]=-0.935;
  F2[1][2][4][4]=-0.612;
  F2[1][2][6][4]=0.655;
  F2[1][2][4][6]=-0.949;
  F2[1][2][6][6]=-0.433;
  F2[1][2][8][6]=0.722;
  F2[1][2][6][8]=-0.940;
  F2[1][2][8][8]=-0.335;
  F2[1][2][10][8]=0.757;
  F2[1][2][8][10]=-0.931;
  F2[1][2][10][10]=-0.274;
  F2[1][2][12][10]=0.778;
  F2[1][2][10][12]=-0.923;
  F2[1][2][12][12]=-0.231;
  F2[1][2][14][12]=0.793;
  F2[1][2][12][14]=-0.917;
  F2[1][2][14][14]=-0.200;
  F2[1][2][16][14]=0.803;
  F2[1][2][14][16]=-0.912;
  F2[1][2][16][16]=-0.177;
  F2[1][2][18][16]=0.811;
  F2[1][2][16][18]=-0.907;
  F2[1][2][18][18]=-0.158;
  F2[1][2][20][18]=0.817;
  F2[1][2][18][20]=-0.904;
  F2[1][2][20][20]=-0.143;
  F2[1][2][22][20]=0.822;

  F2[2][2][2][2]=-0.354;
  F2[2][2][4][2]=0.354;
  F2[2][2][6][2]=-0.101;
  F2[2][2][0][4]=-0.598;
  F2[2][2][2][4]=-0.299;
  F2[2][2][4][4]=0.128;
  F2[2][2][6][4]=0.341;
  F2[2][2][8][4]=-0.171;
  F2[2][2][2][6]=-0.495;
  F2[2][2][4][6]=-0.124;
  F2[2][2][6][6]=0.227;
  F2[2][2][8][6]=0.309;
  F2[2][2][10][6]=-0.206;
  F2[2][2][4][8]=-0.448;
  F2[2][2][6][8]=-0.045;
  F2[2][2][8][8]=0.265;
  F2[2][2][10][8]=0.285;
  F2[2][2][12][8]=-0.228;
  F2[2][2][6][10]=-0.421;
  F2[2][2][8][10]=0.0;
  F2[2][2][10][10]=0.283;
  F2[2][2][12][10]=0.267;
  F2[2][2][14][10]=-0.243;
  F2[2][2][8][12]=-0.403;
  F2[2][2][10][12]=0.029;
  F2[2][2][12][12]=0.294;
  F2[2][2][14][12]=0.253;
  F2[2][2][16][12]=-0.253;
  F2[2][2][10][14]=-0.391;
  F2[2][2][12][14]=0.049;
  F2[2][2][14][14]=0.300;
  F2[2][2][16][14]=0.243;
  F2[2][2][18][14]=-0.261;
  F2[2][2][12][16]=-0.381;
  F2[2][2][14][16]=0.064;
  F2[2][2][16][16]=0.304;
  F2[2][2][18][16]=0.234;
  F2[2][2][20][16]=-0.268;
  F2[2][2][14][18]=-0.374;
  F2[2][2][16][18]=0.075;
  F2[2][2][18][18]=0.307;
  F2[2][2][20][18]=0.227;
  F2[2][2][22][18]=-0.273;
  F2[2][2][16][20]=-0.369;
  F2[2][2][18][20]=0.084;
  F2[2][2][20][20]=0.310;
  F2[2][2][22][20]=0.221;
  F2[2][2][24][20]=-0.277;

  F4[0][4]=-1.069;
  F4[2][4]=0.713;
  F4[4][4]=-0.305;
  F4[6][4]=0.076;
  F4[8][4]=-0.008;
  F4[2][6]=-0.447;
  F4[4][6]=0.670;
  F4[6][6]=-0.447;
  F4[8][6]=0.149;
  F4[10][6]=-0.020;
  F4[4][8]=-0.304;
  F4[6][8]=0.609;
  F4[8][8]=-0.498;
  F4[10][8]=0.194;
  F4[12][8]=-0.030;
  F4[6][10]=-0.243;
  F4[8][10]=0.567;
  F4[10][10]=-0.523;
  F4[12][10]=0.224;
  F4[14][10]=-0.037;
  F4[8][12]=-0.209;
  F4[10][12]=0.537;
  F4[12][12]=-0.537;
  F4[14][12]=0.246;
  F4[16][12]=-0.043;
  F4[10][14]=-0.187;
  F4[12][14]=0.515;
  F4[14][14]=-0.546;
  F4[16][14]=0.263;
  F4[18][14]=-0.048;
  F4[12][16]=-0.173;
  F4[14][16]=0.499;
  F4[16][16]=-0.551;
  F4[18][16]=0.276;
  F4[20][16]=-0.059;
  F4[14][18]=-0.162;
  F4[16][18]=0.486;
  F4[18][18]=-0.555;
  F4[20][18]=0.286;
  F4[22][18]=-0.056;
  F4[16][20]=-0.154;
  F4[18][20]=0.476;
  F4[20][20]=-0.558;
  F4[22][20]=0.295;
  F4[24][20]=-0.059;

// !     Odd spins

  F2[1][1][1][3]=0.500;
  F2[1][1][3][3]=-0.400;
  F2[1][1][5][3]=0.100;
  F2[1][1][3][5]=0.374;
  F2[1][1][5][5]=-0.428;
  F2[1][1][7][5]=0.134;
  F2[1][1][5][7]=0.327;
  F2[1][1][7][7]=-0.436;
  F2[1][1][9][7]=0.153;
  F2[1][1][7][9]=0.303;
  F2[1][1][9][9]=-0.440;
  F2[1][1][11][9]=0.165;
  F2[1][1][9][11]=0.288;
  F2[1][1][11][11]=-0.442;
  F2[1][1][13][11]=0.174;
  F2[1][1][11][13]=0.277;
  F2[1][1][13][13]=-0.444;
  F2[1][1][15][13]=0.180;
  F2[1][1][13][15]=0.270;
  F2[1][1][15][15]=-0.445;
  F2[1][1][17][15]=0.185;
  F2[1][1][15][17]=0.264;
  F2[1][1][17][17]=-0.445;
  F2[1][1][19][17]=0.189;
  F2[1][1][17][19]=0.260;
  F2[1][1][19][19]=-0.446;
  F2[1][1][21][19]=0.192;

  F2[1][2][1][3]=-0.866;
  F2[1][2][3][3]=-0.755;
  F2[1][2][5][3]=0.592;
  F2[1][2][3][5]=-0.949;
  F2[1][2][5][5]=-0.507;
  F2[1][2][7][5]=0.694;
  F2[1][2][5][7]=-0.945;
  F2[1][2][7][7]=-0.378;
  F2[1][2][9][7]=0.742;
  F2[1][2][7][9]=-0.935;
  F2[1][2][9][9]=-0.302;
  F2[1][2][11][9]=0.769;
  F2[1][2][9][11]=-0.927;
  F2[1][2][11][11]=-0.251;
  F2[1][2][13][11]=0.786;
  F2[1][2][11][13]=-0.920;
  F2[1][2][13][13]=-0.215;
  F2[1][2][15][13]=0.798;
  F2[1][2][13][15]=-0.914;
  F2[1][2][15][15]=-0.188;
  F2[1][2][17][15]=0.807;
  F2[1][2][15][17]=-0.910;
  F2[1][2][17][17]=-0.167;
  F2[1][2][19][17]=0.814;
  F2[1][2][17][19]=-0.906;
  F2[1][2][19][19]=-0.150;
  F2[1][2][21][19]=0.820;

  F2[2][2][1][3]=-0.500;
  F2[2][2][5][3]=0.357;
  F2[2][2][7][3]=-0.143;
  F2[2][2][1][5]=-0.535;
  F2[2][2][3][5]=-0.191;
  F2[2][2][5][5]=0.191;
  F2[2][2][7][5]=0.325;
  F2[2][2][9][5]=-0.191;
  F2[2][2][3][7]=-0.468;
  F2[2][2][5][7]=-0.078;
  F2[2][2][7][7]=0.249;
  F2[2][2][9][7]=0.296;
  F2[2][2][11][7]=-0.218;
  F2[2][2][5][9]=-0.433;
  F2[2][2][7][9]=-0.020;
  F2[2][2][9][9]=0.275;
  F2[2][2][11][9]=0.275;
  F2[2][2][13][9]=-0.236;
  F2[2][2][7][11]=-0.411;
  F2[2][2][9][11]=0.016;
  F2[2][2][11][11]=0.289;
  F2[2][2][13][11]=0.260;
  F2[2][2][15][11]=-0.248;
  F2[2][2][9][13]=-0.396;
  F2[2][2][11][13]=0.040;
  F2[2][2][13][13]=0.297;
  F2[2][2][15][13]=0.248;
  F2[2][2][17][13]=-0.258;
  F2[2][2][11][15]=-0.386;
  F2[2][2][13][15]=0.057;
  F2[2][2][15][15]=0.302;
  F2[2][2][17][15]=0.238;
  F2[2][2][19][15]=-0.265;
  F2[2][2][13][17]=-0.378;
  F2[2][2][15][17]=0.070;
  F2[2][2][17][17]=0.306;
  F2[2][2][19][17]=0.231;
  F2[2][2][21][17]=-0.270;
  F2[2][2][15][19]=-0.371;
  F2[2][2][17][19]=0.080;
  F2[2][2][19][19]=0.309;
  F2[2][2][21][19]=0.224;
  F2[2][2][23][19]=-0.275;

  F4[1][5]=-0.617;
  F4[3][5]=0.705;
  F4[5][5]=-0.397;
  F4[7][5]=0.118;
  F4[9][5]=-0.015;
  F4[3][7]=-0.358;
  F4[5][7]=0.637;
  F4[7][7]=-0.478;
  F4[9][7]=0.174;
  F4[11][7]=-0.025;
  F4[5][9]=-0.268;
  F4[7][9]=0.586;
  F4[9][9]=-0.512;
  F4[11][9]=0.210;
  F4[13][9]=-0.034;
  F4[7][11]=-0.224;
  F4[9][11]=0.551;
  F4[11][11]=-0.531;
  F4[13][11]=0.236;
  F4[15][11]=-0.041;
  F4[9][13]=-0.197;
  F4[11][13]=0.525;
  F4[13][13]=-0.542;
  F4[15][13]=0.255;
  F4[17][13]=-0.046;
  F4[11][15]=-0.179;
  F4[13][15]=0.507;
  F4[15][15]=-0.549;
  F4[17][15]=0.270;
  F4[19][15]=-0.051;
  F4[13][17]=-0.167;
  F4[15][17]=0.492;
  F4[17][17]=-0.554;
  F4[19][17]=0.281;
  F4[21][17]=-0.054;
  F4[15][19]=-0.158;
  F4[17][19]=0.481;
  F4[19][19]=-0.557;
  F4[21][19]=0.291;
  F4[23][19]=-0.058;

  if (rank == 2 ) {
    return F2[L1][L2][(int)(2*spin1)][(int)(2*spin2)];
  }
  else if ((rank == 4)&&(L1==2)&&(L2==2)){
    return F4[(int)(2*spin1)][(int)(2*spin2)];
  }
  else{
    return 0.0;
  }
}
