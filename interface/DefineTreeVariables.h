//This is the list and types of the variables saved in the TTree;
//New variables must be declared here

#include "TString.h"

//event variables
int runnumber;
int lumiblock;
int eventNumber;
float timestamp;
int bunchXing;
int orbitNum;

//primary vertex
float PV_x;
float PV_y;
float PV_z;

float PV_xxE;
float PV_yyE;
float PV_zzE;
float PV_xyE;
float PV_xzE;
float PV_yzE;

float PV_normchi2;
float PV_Nvtx;

//luminosity
float lumiperblock;
float beam1Intensity;
float beam2Intensity;

// HLT
std::vector<TString> hlt_path;

//digi variables
std::vector<short> digi_wheel;
std::vector<short> digi_sector;
std::vector<short> digi_station;
std::vector<short> digi_sl;
std::vector<short> digi_layer;
std::vector<short> digi_wire;
std::vector<float> digi_time;

//DT segment variables
std::vector<short> segm4D_wheel;
std::vector<short> segm4D_sector;
std::vector<short> segm4D_station;

std::vector<short> segm4D_hasPhi;
std::vector<short> segm4D_hasZed;

std::vector<float> segm4D_x_pos_loc;
std::vector<float> segm4D_y_pos_loc;
std::vector<float> segm4D_z_pos_loc;
std::vector<float> segm4D_x_dir_loc;
std::vector<float> segm4D_y_dir_loc;
std::vector<float> segm4D_z_dir_loc;

std::vector<float> segm4D_cosx;
std::vector<float> segm4D_cosy;
std::vector<float> segm4D_cosz;
std::vector<float> segm4D_phi;
std::vector<float> segm4D_theta;
std::vector<float> segm4D_eta;

std::vector<float> segm4D_t0;
std::vector<float> segm4D_vdrift;
std::vector<float> segm4D_phinormchi2;
std::vector<short> segm4D_phinhits;

std::vector<float> segm4D_znormchi2;
std::vector<short> segm4D_znhits;

TClonesArray *segm4D_phiHits_Pos;
TClonesArray *segm4D_phiHits_PosCh;
TClonesArray *segm4D_phiHits_PosErr;
TClonesArray *segm4D_phiHits_Side;
TClonesArray *segm4D_phiHits_Wire;
TClonesArray *segm4D_phiHits_Layer;
TClonesArray *segm4D_phiHits_SuperLayer;
TClonesArray *segm4D_phiHits_Time;
TClonesArray *segm4D_phiHits_TimeCali;

TClonesArray *segm4D_hitsExpPos;
TClonesArray *segm4D_hitsExpWire;

TClonesArray *segm4D_zHits_Pos;
TClonesArray *segm4D_zHits_PosCh;
TClonesArray *segm4D_zHits_PosErr;
TClonesArray *segm4D_zHits_Side;
TClonesArray *segm4D_zHits_Wire;
TClonesArray *segm4D_zHits_Layer;
TClonesArray *segm4D_zHits_Time;
TClonesArray *segm4D_zHits_TimeCali;

//CSC segment variables
std::vector<short> cscsegm_ring;
std::vector<short> cscsegm_chamber;
std::vector<short> cscsegm_station;
std::vector<float> cscsegm_cosx;
std::vector<float> cscsegm_cosy;
std::vector<float> cscsegm_cosz;
std::vector<float> cscsegm_phi;
std::vector<float> cscsegm_eta;
std::vector<float> cscsegm_normchi2;
std::vector<short> cscsegm_nRecHits;

//TM variables
std::vector<short> ltTwinMuxIn_wheel;
std::vector<short> ltTwinMuxIn_sector;
std::vector<short> ltTwinMuxIn_station;
std::vector<short> ltTwinMuxIn_quality;
std::vector<short> ltTwinMuxIn_bx;
std::vector<float> ltTwinMuxIn_phi;
std::vector<float> ltTwinMuxIn_phiB;
std::vector<short> ltTwinMuxIn_is2nd;

std::vector<short> ltTwinMuxOut_wheel;
std::vector<short> ltTwinMuxOut_sector;
std::vector<short> ltTwinMuxOut_station;
std::vector<short> ltTwinMuxOut_quality;
std::vector<short> ltTwinMuxOut_rpcbit;
std::vector<short> ltTwinMuxOut_bx;
std::vector<float> ltTwinMuxOut_phi;
std::vector<float> ltTwinMuxOut_phiB;
std::vector<short> ltTwinMuxOut_is2nd;

std::vector<short> ltTwinMux_thBx;
std::vector<short> ltTwinMux_thWheel;
std::vector<short> ltTwinMux_thSector;
std::vector<short> ltTwinMux_thStation;
std::vector<short> ltTwinMux_thHits;

//muon variables
std::vector<short> STAMu_isMuGlobal;
std::vector<short> STAMu_isMuTracker;
std::vector<int>   STAMu_numberOfChambers;
std::vector<int>   STAMu_numberOfMatches;
std::vector<int>   STAMu_numberOfHits;
std::vector<int>   STAMu_segmIndex;

std::vector<float> Mu_px_mu;
std::vector<float> Mu_py_mu;
std::vector<float> Mu_pz_mu;
std::vector<float> Mu_phi_mu;
std::vector<float> Mu_eta_mu;
std::vector<short> STAMu_recHitsSize;
std::vector<float> STAMu_normchi2Mu;
std::vector<short> STAMu_chargeMu;
std::vector<float> STAMu_dxyMu;
std::vector<float> STAMu_dzMu;

std::vector<float> GLBMu_normchi2Mu;
std::vector<float> GLBMu_dxyMu;
std::vector<float> GLBMu_dzMu;

std::vector<int> GLBMu_numberOfPixelHits;
std::vector<int> GLBMu_numberOfTrackerHits;

std::vector<float> GLBMu_tkIsoR03;
std::vector<float> GLBMu_ntkIsoR03;
std::vector<float> GLBMu_emIsoR03;
std::vector<float> GLBMu_hadIsoR03;

std::vector<float> STAMu_caloCompatibility;

std::vector<float> STAMu_z_mb2;
std::vector<float> STAMu_phi_mb2;
std::vector<float> STAMu_pseta_mb2;

//GMT
std::vector<short> gmt_bx;
std::vector<float> gmt_phi;
std::vector<float> gmt_eta;
std::vector<float> gmt_pt;
std::vector<short> gmt_qual;
std::vector<short> gmt_charge;
std::vector<int>   gmt_tf_idx;

//GT // legacy
std::vector<short> gt_algo_bx;
std::vector<short> gt_algo_bit;
std::vector<short> gt_tt_bx;
std::vector<short> gt_tt_bit;

//RPC
std::vector<int>   rpc_region;
std::vector<int>   rpc_clusterSize;
std::vector<int>   rpc_strip;
std::vector<int>   rpc_bx;
std::vector<int>   rpc_station;
std::vector<int>   rpc_sector;
std::vector<int>   rpc_layer;
std::vector<int>   rpc_subsector;
std::vector<int>   rpc_roll;
std::vector<int>   rpc_ring;

int Bmtf_Size;
std::vector<short int> Bmtf_Pt;
std::vector<short int> Bmtf_Eta;
std::vector<short int> Bmtf_Phi;
std::vector<short int> Bmtf_GlobalPhi;
std::vector<short int> Bmtf_qual;
std::vector<short int> Bmtf_ch;
std::vector<short int> Bmtf_bx;
std::vector<short int> Bmtf_processor;
std::vector<short int> Bmtf_trAddress;
std::vector<short int> Bmtf_wh;
std::vector<short int> Bmtf_FineBit;

int Bmtf_phSize;
std::vector<int> Bmtf_phBx;
std::vector<int> Bmtf_phWh;
std::vector<int> Bmtf_phSe;
std::vector<int> Bmtf_phSt;
std::vector<float>  Bmtf_phAng;
std::vector<float>  Bmtf_phBandAng;
std::vector<int> Bmtf_phCode;
std::vector<int> Bmtf_phTs2Tag;

int Bmtf_thSize;
std::vector<int>   Bmtf_thBx;
std::vector<int>   Bmtf_thWh;
std::vector<int>   Bmtf_thSe;
std::vector<int>   Bmtf_thSt;
std::vector<int> Bmtf_thTheta;
std::vector<int> Bmtf_thCode; 

std::vector<int> RpcDigi_TwinMux_bx;
std::vector<int> RpcDigi_TwinMux_strip;
std::vector<int> RpcDigi_TwinMux_region;
std::vector<int> RpcDigi_TwinMux_ring;
std::vector<int> RpcDigi_TwinMux_station;
std::vector<int> RpcDigi_TwinMux_layer;
std::vector<int> RpcDigi_TwinMux_sector;
std::vector<int> RpcDigi_TwinMux_subsector;
std::vector<int> RpcDigi_TwinMux_roll;
std::vector<int> RpcDigi_TwinMux_trIndex;
std::vector<int> RpcDigi_TwinMux_det;
std::vector<int> RpcDigi_TwinMux_subdetId;
std::vector<int> RpcDigi_TwinMux_rawId;

std::vector<int> RpcRechit_TwinMux_region;
std::vector<int> RpcRechit_TwinMux_clusterSize;
std::vector<int> RpcRechit_TwinMux_strip;
std::vector<int> RpcRechit_TwinMux_bx;
std::vector<int> RpcRechit_TwinMux_station;
std::vector<int> RpcRechit_TwinMux_sector;
std::vector<int> RpcRechit_TwinMux_layer;
std::vector<int> RpcRechit_TwinMux_subsector;
std::vector<int> RpcRechit_TwinMux_roll;
std::vector<int> RpcRechit_TwinMux_ring;
std::vector<float> RpcRechit_TwinMux_Loc_x;
std::vector<float> RpcRechit_TwinMux_Loc_y;
std::vector<float> RpcRechit_TwinMux_Loc_z;
std::vector<float> RpcRechit_TwinMux_Glob_x;
std::vector<float> RpcRechit_TwinMux_Glob_y;
std::vector<float> RpcRechit_TwinMux_Glob_z;
std::vector<float> RpcRechit_TwinMux_Glob_eta;
std::vector<float> RpcRechit_TwinMux_Glob_phi;


std::vector<float> DT_extrapolated_OnRPC_BX;
std::vector<float> DT_extrapolated_OnRPC_Loc_x;
std::vector<float> DT_extrapolated_OnRPC_Loc_y;
std::vector<float> DT_extrapolated_OnRPC_Loc_z;
std::vector<float> DT_extrapolated_OnRPC_Glob_x;
std::vector<float> DT_extrapolated_OnRPC_Glob_y;
std::vector<float> DT_extrapolated_OnRPC_Glob_z;
std::vector<float> DT_extrapolated_OnRPC_Glob_eta;
std::vector<float> DT_extrapolated_OnRPC_Glob_phi;
std::vector<float> DT_extrapolated_OnRPC_Region;
std::vector<float> DT_extrapolated_OnRPC_Sector;
std::vector<float> DT_extrapolated_OnRPC_Station;
std::vector<float> DT_extrapolated_OnRPC_Layer;
std::vector<float> DT_extrapolated_OnRPC_Roll;
std::vector<float> DT_extrapolated_OnRPC_Ring;
std::vector<float> DT_extrapolated_OnRPC_Stripw;
