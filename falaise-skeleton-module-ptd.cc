#include <bayeux/dpp/chain_module.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/tracker_trajectory_solution.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/unified_digitized_data.h>
#include <falaise/snemo/datamodels/tracker_cluster.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/geomid_utils.h>  
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>
#include <bayeux/geomtools/geom_info.h>
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/locator_plugin.h>
#include <falaise/snemo/geometry/xcalo_locator.h>
#include <geomtools/geometry_service.h>
#include <datatools/service_manager.h>
#include <geomtools/manager.h>
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include "TLatex.h"
#include "TVector3.h"
#include "TParameter.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <cmath>

using namespace std;


class falaise_skeleton_module_ptd : public dpp::chain_module
{

public:
  // Constructor
  falaise_skeleton_module_ptd();

  // Destructor
  virtual ~falaise_skeleton_module_ptd();

  // Initialisation function
  virtual void initialize (const datatools::properties &,
                           datatools::service_manager & serviceManager,
			   dpp::module_handle_dict_type &);
  
  double compute_ellipse(double y_vertex, double z_vertex, int& source_num, double &dy, double &dz) const;
  double extract_unix_start_time(int run_number) const;
  // Event processing function
  dpp::chain_module::process_status process (datatools::things & event);
  
private:
  int event_number;
  int  ptd_event_counter, nb_unfitted_cluster, nb_unfitted_cells;
  int cellules_non_associated, cellules_SD_non_associated, number_of_kinks, number_of_electrons;
  bool ptd_details, has_an_electron_and_positron, has_SD_electron_and_positron, has_SD_two_electrons, opposite_side_e_gamma, same_side_elec, same_cluster_elec, has_a_same, ttd_details, sd_calo_details, sd_tracker_details, hit_the_same_calo_hit;

  double diff_time_elec, angle_3D_between_ep_em, delta_y_elec, delta_z_elec, energy_elec_sum, corrected_energy_elec_sum, time_of_flight_gamma, internal_theoretical_time_diff, external_theoretical_time_diff,t1_th, t2_th, start_run_time, end_run_time, delta_r_calo, kink_angle, one_kink_x, one_kink_y, one_kink_z, unix_start_time=0, first_time, last_time, clustering_distance;
  int nb_gamma, nb_elec_ptd_per_event,nb_elec_SD_per_event, nb_gamma_per_elec, nb_elec_linked;
  vector<string>  type_elec, g4_process, material, vertex_type;
  vector<int> gamma_type, num_om, track_number, num_source, num_om_gamma, num_om_gamma_f, nb_gamma_per_elec_vec, nb_elec_linked_vec;
  vector<double> time_gamma, time_gamma_before, time_gamma_after, angle_SD, charge_elec, ellipse_source, delta_y_source, delta_z_source, charge_gamma, total_energy_gamma, total_energy_gamma_vis, total_energy_electron;
  vector<double> energy_gamma, energy_elec, corrected_energy_elec, time_elec,track_lenght, energy_gamma_after, vertex_SD_x, vertex_SD_y, vertex_SD_z , delta_y_calo_center, delta_z_calo_center, evis_bcu, evis_bc, evis_u, evis, chi2_track, edep_bcu, gamma_evis_bcu, gamma_evis_bc, gamma_evis_u, gamma_evis, gamma_edep, gamma_edep_bcu, gamma_evis_bcu_vec, gamma_evis_bc_vec, gamma_evis_u_vec, gamma_evis_vec, gamma_edep_vec, gamma_edep_bcu_vec;
  vector<double> vertex_3D_start_x, vertex_3D_start_y, vertex_3D_start_z,vertex_3D_end_x, vertex_3D_end_y, vertex_3D_end_z, vertex_gamma, kink_x, kink_y, kink_z, gamma_vertex_x, gamma_vertex_y, gamma_vertex_z, dt_gamma_elec, dt_gamma_elec_c, base_dt_gamma_elec, charge_gamma_finale, charge_elec_finale;
  vector<int> side_elec, number_of_kinks_per_track;
  const geomtools::manager* _geo_manager_ = nullptr;
  std::map<int, double> y_map;
  std::map<int, double> z_map;
  std::map<int,double> calib_param;
  TFile *save_file;
  TTree *tree;
  //data source positions
  std::vector<double> y_source_pos = {
    -2087.5, -2087.5, -2087.5, -2087.5, -2087.5, -2087.5, -2087.5,
    -1252.5, -1252.5, -1252.5, -1252.5, -1252.5, -1252.5, -1252.5,
    -417.5,  -417.5,  -417.5,  -417.5,  -417.5,  -417.5,  -417.5,
    417.5,   417.5,   417.5,   417.5,   417.5,   417.5,   417.5,
    1252.5,  1252.5,  1252.5,  1252.5,  1252.5,  1252.5,  1252.5,
    2087.5,  2087.5,  2087.5,  2087.5,  2087.5,  2087.5,  2087.5
  };
  
  std::vector<double> z_source_pos = {
    1317.5,   882.5,   447.5,     2.5,  -442.5,  -882.5, -1317.5,
    1317.5,   887.5,   447.5,     2.5,  -442.5,  -887.5, -1317.5,
    1312.5,   882.5,   447.5,     2.5,  -442.5,  -882.5, -1317.5,
    1317.5,   887.5,   447.5,     2.5,  -437.5,  -882.5, -1317.5,
    1317.5,   882.5,   442.5,     2.5,  -442.5,  -882.5, -1317.5,
    1317.5,   882.5,   442.5,     2.5,  -442.5,  -882.5, -1317.5
  };

  // //simu source position
  // std::vector<double> y_source_pos = {
  //   -2075, -2075, -2075, -2075, -2075, -2075, -2075,
  //   -1225, -1225, -1225, -1225, -1225, -1225, -1225,
  //   -425,  -425,  -425,  -425,  -425,  -425,  -425,
  //   425,   425,   425,   425,   425,   425,   425,
  //   1275,  1275,  1275,  1275,  1275,  1275,  1275,
  //   2075,  2075,  2075,  2075,  2075,  2075,  2075
  // };
  
  // std::vector<double> z_source_pos = {
  //   1305,   885,   425,     15,  -405,  -855, -1275,
  //   1305,   885,   425,     15,  -405,  -855, -1275,
  //   1305,   885,   425,     15,  -405,  -855, -1275,
  //   1305,   885,   425,     15,  -405,  -855, -1275,
  //   1305,   885,   425,     15,  -405,  -855, -1275,
  //   1305,   885,   425,     15,  -405,  -855, -1275
  // };

  

  std::vector<double> source_pos_num = {
    6,  5,  4,  3,   2,  1,  0,
    13, 12, 11, 10,  9,  8,  7,
    20, 19, 18, 17, 16, 15, 14,
    27, 26, 25, 24, 23, 22, 21,
    34, 33, 32, 31, 30, 29, 28,
    41, 40, 39, 38, 37, 36, 35,    
  };
// std::vector<double> source_pos_num = {
//   9, 16, 23, 7, 8, 15, 30,
//   2, 21, 22, 29, 14, 28, 37,
//   1, 0, 36, 35
// };

// std::vector<double> y_source_pos = {
//   -1255.48, -416.985, 415.925, -1255.48, -1255.48, -416.985, 1254.43,
//   -2088.39, 415.925, 415.925, 1254.43, -416.985, 1254.43, 2087.34,
//   -2088.39, -2088.39, 2087.34, 2087.34
// };

// std::vector<double> z_source_pos = {
//   1139.53, 1132.11, 1135.82, 260.265, 701.755, 698.045, 1132.11,
//   1135.82, 256.555, 698.045, 694.335, 260.265, 252.845, 1132.11,
//   698.045, 260.265, 698.045, 252.845
// }; 
  
  // Macro to register the module
  DPP_MODULE_REGISTRATION_INTERFACE(falaise_skeleton_module_ptd);

};

////////////////////////////////////////////////////////////////////

// Macro to add the module in the global register of data processing modules:
// The module defined by this class 'falaise_skeleton_module_ptd' will be registered
// with the label ID 'FalaiseSkeletonModule_PTD' (to use in pipeline configuration file)
DPP_MODULE_REGISTRATION_IMPLEMENT(falaise_skeleton_module_ptd, "FalaiseSkeletonModule_PTD")


falaise_skeleton_module_ptd::falaise_skeleton_module_ptd()
{
  save_file = new TFile("extracted_data.root", "RECREATE");
  tree = new TTree("Event", "Event information");  
  tree->Branch("event_number", &event_number);
  tree->Branch("nb_gamma", &nb_gamma);
  tree->Branch("vertex_3D_start_x", &vertex_3D_start_x);
  tree->Branch("vertex_3D_start_y", &vertex_3D_start_y);
  tree->Branch("vertex_3D_start_z", &vertex_3D_start_z);
  tree->Branch("ellipse_source", &ellipse_source);
  tree->Branch("delta_y_source", &delta_y_source);
  tree->Branch("delta_z_source", &delta_z_source);
  tree->Branch("num_source", &num_source);
  tree->Branch("vertex_3D_end_x", &vertex_3D_end_x);
  tree->Branch("vertex_3D_end_y", &vertex_3D_end_y);
  tree->Branch("vertex_3D_end_z", &vertex_3D_end_z);
  tree->Branch("delta_y_calo_center", &delta_y_calo_center);
  tree->Branch("delta_z_calo_center", &delta_z_calo_center);
  tree->Branch("num_om",&num_om);
  tree->Branch("num_om_gamma",&num_om_gamma);
  tree->Branch("num_om_gamma_f",&num_om_gamma_f);
  tree->Branch("charge_gamma",&charge_gamma);
  tree->Branch("energy_gamma",&energy_gamma);
  tree->Branch("total_energy_gamma",&total_energy_gamma);
  tree->Branch("total_energy_gamma_vis",&total_energy_gamma_vis);
  tree->Branch("total_energy_electron",&total_energy_electron);
  tree->Branch("time_gamma",&time_gamma);
  tree->Branch("gamma_type",&gamma_type);
  tree->Branch("gamma_vertex_x",&gamma_vertex_x);
  tree->Branch("gamma_vertex_y",&gamma_vertex_y);
  tree->Branch("gamma_vertex_z",&gamma_vertex_z);
  tree->Branch("nb_gamma_per_elec",&nb_gamma_per_elec_vec);
  tree->Branch("nb_elec_linked",&nb_elec_linked_vec);
  tree->Branch("charge_elec_finale",&charge_elec_finale);
  tree->Branch("charge_gamma_finale",&charge_gamma_finale);
  tree->Branch("dt_gamma_elec_real_E",&dt_gamma_elec);
  tree->Branch("dt_gamma_elec_c",&dt_gamma_elec_c);
  tree->Branch("base_dt_gamma_elec",&base_dt_gamma_elec);
  tree->Branch("nb_elec_ptd_per_event", &nb_elec_ptd_per_event);
  tree->Branch("number_of_kinks_per_track", &number_of_kinks_per_track);
  tree->Branch("energy_elec", &energy_elec);
  tree->Branch("evis", &evis);
  tree->Branch("evis_bcu", &evis_bcu);
  tree->Branch("edep_bcu", &edep_bcu);
  tree->Branch("evis_bc", &evis_bc);
  tree->Branch("evis_u", &evis_u);
  tree->Branch("gamma_evis", &gamma_evis_vec);
  tree->Branch("gamma_evis_bcu", &gamma_evis_bcu_vec);
  tree->Branch("gamma_evis_bc", &gamma_evis_bc_vec);
  tree->Branch("gamma_evis_u", &gamma_evis_u_vec);
  tree->Branch("gamma_edep", &gamma_edep_vec);
  tree->Branch("gamma_edep_bcu", &gamma_edep_bcu_vec);
  tree->Branch("charge_elec", &charge_elec);
  tree->Branch("time_elec", &time_elec);
  tree->Branch("side_elec", &side_elec);
  tree->Branch("type_elec", &type_elec);
  tree->Branch("unix_start_time", &unix_start_time);
  tree->Branch("energy_elec_sum", &energy_elec_sum);
  tree->Branch("cellules_non_associated", &cellules_non_associated);
  tree->Branch("number_of_electrons",&number_of_electrons);
  tree->Branch("clustering_distance",&clustering_distance);
  tree->Branch("nb_unfitted_cluster",&nb_unfitted_cluster);
  tree->Branch("nb_unfitted_cells",&nb_unfitted_cells);
  tree->Branch("chi2_track",&chi2_track);

  std::ifstream file_calo("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/OM_positions.csv");  
  int om_num;
  double y, z;
  
  while (file_calo >> om_num >> y >> z) {
    y_map[om_num] = y;
    z_map[om_num] = z;
  }
}


falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd()
{
  double time = last_time - first_time; 
  save_file->cd();
  TParameter<double> param("run_time", time);
  param.Write();
  tree->Write();
  save_file->Close();
  delete save_file;
  std::cout << "falaise_skeleton_module_ptd::~falaise_skeleton_module_ptd() called" << std::endl;
}


//read calibration param
std::map<int,double> read_calib_param()
{
    std::map<int,double> om_map;

    std::ifstream file("/sps/nemo/scratch/granjon/full_gain_analysis/Bi/root/compute_a/calibration_phase_2b.txt");
    std::string line;
    bool header_skipped = false;
    while (std::getline(file, line))
    {
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }
        std::istringstream iss(line);
        int om;
        double value;
        if (!(iss >> om >> value))
            continue; 
        om_map[om] = value;
    }

    file.close();
    return om_map;
}



void falaise_skeleton_module_ptd::initialize (const datatools::properties & module_properties, datatools::service_manager & serviceManager, dpp::module_handle_dict_type &)
{
  std::cout << "falaise_skeleton_module_ptd::initialize() called" << std::endl;
  event_number=0;
  first_time=0;
  last_time=0;

  if ( module_properties.has_key("ptd_details"))
    ptd_details = module_properties.fetch_boolean("ptd_details");
  else ptd_details = false;

  if (module_properties.has_key("sigma_r")) {
    module_properties.fetch("sigma_r", clustering_distance);
  }

  if ( module_properties.has_key("calo_details"))
    sd_calo_details = module_properties.fetch_boolean("calo_details");
  else sd_calo_details = false;

  if ( module_properties.has_key("tracker_details"))
    sd_tracker_details = module_properties.fetch_boolean("tracker_details");
  else sd_tracker_details = false;

  _geo_manager_ = &serviceManager.get<geomtools::geometry_service>("geometry").get_geom_manager();

  calib_param = read_calib_param();

  this->_set_initialized(true);
}

double snemo_run_time (int run_number)
{
  const char *cbd_base_path = "/sps/nemo/snemo/snemo_data/raw_data/CBD";
  std::vector<std::string> log_paths;
  log_paths.push_back(Form("%s/run-%d/snemo_trigger_run-%d.log", cbd_base_path, run_number, run_number));
  for (int crate=6; crate>=0; crate--)
    log_paths.push_back(Form("%s/run-%d/snemo_crate-%d_run-%d.log", cbd_base_path, run_number, crate, run_number));
  int unixtime;
  int run_start=0;
  for (const std::string & log_path : log_paths)
    {
      std::ifstream log_file (log_path);
      if (!log_file.is_open()) continue;
      std::string log_line;
      while (getline(log_file, log_line))
        {
          size_t unixtime_index = log_line.find("run.run_unixtime_ms=");
          if (unixtime_index == std::string::npos)
            continue;
          unixtime = std::stoi(log_line.substr(20));
          if (unixtime > run_start)
            run_start = unixtime;
	    }
    }
  return unixtime;
}





double falaise_skeleton_module_ptd::extract_unix_start_time(int run_number) const
{
    std::string cmd = "grep ^" + std::to_string(run_number) + " /sps/nemo/scratch/chauveau/commissioning/software/run-sync-time.txt | awk '{print $2}'";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
	std::cerr << "Erreur: impossible d'exécuter grep\n";
	return -1;
    }

    double timestamp = -1;
    fscanf(pipe, "%lf", &timestamp);
    pclose(pipe);
    return timestamp;
}


double falaise_skeleton_module_ptd::compute_ellipse(double y_vertex, double z_vertex, int& source_num, double &best_dy, double &best_dz) const
{
  best_dy = 0.0;
  best_dz = 0.0;
  
  double min_dist2 = std::numeric_limits<double>::max();

  for (size_t i = 0; i < y_source_pos.size(); ++i) {
    double dy = y_source_pos[i] - y_vertex;
    double dz = z_source_pos[i] - z_vertex;
    double dist2 = (dy * dy) / (25.0 * 25.0) + (dz * dz) / (30.0 * 30.0);

    if (dist2 < min_dist2) {
      min_dist2 = dist2;
      source_num = source_pos_num[i];
      best_dy = dy;
      best_dz = dz;
    }
  }

  return min_dist2; // distance minimale
}




dpp::chain_module::process_status falaise_skeleton_module_ptd::process (datatools::things & event)
{
  // Skip processing if PTD bank is not present
  if (!event.has("PTD"))
    {
      std::cout << "======== no PTD bank in event " << ptd_event_counter++ << " ========" << std::endl;
      return dpp::base_module::PROCESS_SUCCESS;
      event_number++;
    }

  // Retrieve the PTD bank  
  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");

  nb_gamma=0;
  vertex_gamma.clear();
  gamma_vertex_x.clear();
  nb_gamma_per_elec_vec.clear();
  nb_elec_linked_vec.clear();
  gamma_vertex_y.clear();
  gamma_vertex_z.clear();
  number_of_kinks_per_track.clear();
  nb_elec_ptd_per_event=0;
  nb_elec_SD_per_event=0;
  cellules_non_associated=0;
  cellules_SD_non_associated=0;
  gamma_type.clear();
  time_gamma.clear();
  time_gamma_before.clear();
  time_gamma_after.clear();
  energy_gamma_after.clear();
  type_elec.clear();
  side_elec.clear();
  energy_elec.clear();
  num_om.clear();
  dt_gamma_elec.clear();
  dt_gamma_elec_c.clear();
  base_dt_gamma_elec.clear();
  num_om_gamma.clear();
  num_om_gamma_f.clear();
  corrected_energy_elec.clear();
  time_elec.clear();
  track_lenght.clear();
  energy_gamma.clear();
  vertex_3D_start_x.clear();
  vertex_3D_start_y.clear();
  vertex_3D_start_z.clear();
  vertex_3D_end_x.clear();
  edep_bcu.clear();
  evis_bcu.clear();
  evis_bc.clear();
  evis_u.clear();
  evis.clear();
  gamma_evis_bcu.clear();
  gamma_evis_bc.clear();
  gamma_evis_u.clear();
  gamma_evis.clear();
  gamma_edep.clear();
  gamma_edep_bcu.clear();
  gamma_evis_bcu_vec.clear();
  gamma_evis_bc_vec.clear();
  gamma_evis_u_vec.clear();
  gamma_evis_vec.clear();
  gamma_edep_vec.clear();
  gamma_edep_bcu_vec.clear();
  vertex_3D_end_y.clear();
  vertex_3D_end_z.clear();
  ellipse_source.clear();
  delta_y_source.clear();
  delta_z_source.clear();
  vertex_SD_x.clear();
  vertex_SD_y.clear();
  vertex_SD_z.clear();
  kink_x.clear();
  kink_y.clear();
  kink_z.clear();
  time_of_flight_gamma = 0.0;
  energy_elec_sum=0;
  corrected_energy_elec_sum = 0;
  nb_gamma_per_elec=0;
  nb_elec_linked=0;
  kink_angle=0.0;
  same_side_elec=0;
  same_cluster_elec = 0;
  diff_time_elec= 0;
  delta_y_elec=0;
  delta_z_elec=0;
  number_of_electrons=0;
  total_energy_electron.clear();
  total_energy_gamma.clear();
  total_energy_gamma_vis.clear();
  angle_3D_between_ep_em = 0.0;
  angle_SD.clear();
  g4_process.clear();
  material.clear();
  vertex_type.clear();
  charge_elec.clear();
  charge_gamma.clear();
  track_number.clear();
  delta_r_calo=0.0;
  has_an_electron_and_positron = 0;
  has_SD_two_electrons=0;
  has_SD_electron_and_positron=0;
  opposite_side_e_gamma=0;
  internal_theoretical_time_diff = 0.0;
  external_theoretical_time_diff =0.0;
  t1_th=0.0;
  t2_th =0.0;
  start_run_time=0.0;
  end_run_time=0.0;
  one_kink_x=0.0;
  one_kink_y=0.0;
  one_kink_z=0.0;
  nb_unfitted_cluster=0;
  nb_unfitted_cells=0;
  number_of_kinks=0;
  chi2_track.clear();
  num_source.clear();
  delta_y_calo_center.clear();
  delta_z_calo_center.clear();
  charge_gamma_finale.clear();
  charge_elec_finale.clear();
  hit_the_same_calo_hit=false;
  // PTD extraction
  if (ptd_details)    
    {
      const snemo::datamodel::event_header & header = event.get<snemo::datamodel::event_header>("EH"); 
      const datatools::event_id & eh_id = header.get_id();
      const int run_id  = eh_id.get_run_number();
      
      if(header.has_timestamp()){
	double time = (header.get_timestamp().get_seconds()+header.get_timestamp().get_picoseconds()/1E12)*1e-9;
	last_time = time;
	if(first_time==0){
	  first_time = time;
	}
      }
      if(unix_start_time==0){
	if(run_id<1556){
	  unix_start_time = snemo_run_time(run_id);
	}
	else{
	  unix_start_time=extract_unix_start_time(run_id);
	}
      }
    
      const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");
       
      std::vector<std::vector<int>> same_clusters;
      std::vector<std::vector<double>> kinks_vec;
      int nb_kinks_save=0;
      int count_equal_clusters = 0;
      int nb_clusters=0;
      //check unfitted cluster
      
      const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");

      const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();
      for (const auto & trajectory : ttd_solution.get_trajectories()) {
	if(trajectory->get_fit_infos().has_chi2() && trajectory->get_fit_infos().has_ndof()){
	  //cout<<trajectory->get_fit_infos().get_chi2()<<endl;
	  chi2_track.push_back(trajectory->get_fit_infos().get_chi2() / trajectory->get_fit_infos().get_ndof());
	}
	if (ttd_solution.has_unfitted_clusters()){
	  nb_unfitted_cluster=ttd_solution.get_unfitted_clusters().size();
	  cellules_non_associated++;	
	}
      }

      const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");
	const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();
	nb_unfitted_cells=tcd_solution.get_unclustered_hits().size();
	


      
	//check gammas
	vector<int>  particle_type;
	if (PTD.hasIsolatedCalorimeters()) {//gamma start
	  const snemo::datamodel::CalorimeterHitHdlCollection &cc_collection = PTD.isolatedCalorimeters();
	  for (const auto &it_hit : cc_collection) {
	    //possible de tout récuperer mais pas de vertex -> pas de gammatype
	    const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & pcd_calo_hit = pCD.calorimeter_hits()[it_hit->get_hit_id()];
	    //energy_gamma.push_back(it_hit->get_energy());
	    charge_gamma.push_back((-pcd_calo_hit->get_charge())/(1E-9*CLHEP::volt*CLHEP::second * 0.00023841858));            
	    int om_number = snemo::datamodel::om_num(it_hit->get_geom_id());
	    const geomtools::manager & geoMgr = *_geo_manager_;
	    //const geomtools::manager & geoMgr = serviceManager.get<geomtools::geometry_service>("geometry").get_geom_manager();
	    
	    uint32_t moduleNumber = 0;

	    int type = it_hit->get_geom_id().get_type();
	    geomtools::vector_3d pos;
	    if(type == 1302){
	      snemo::geometry::calo_locator CL(moduleNumber, geoMgr,{});
	      pos = CL.getBlockPosition(it_hit->get_geom_id());
	    }
	    else if(type == 1232){
	      snemo::geometry::xcalo_locator CL(moduleNumber, geoMgr,{});
	      pos = CL.getBlockPosition(it_hit->get_geom_id());
	    }
	    else if(type == 1251){
	      snemo::geometry::gveto_locator CL(moduleNumber, geoMgr,{});
	      pos = CL.getBlockPosition(it_hit->get_geom_id());
	    }
	    gamma_vertex_x.push_back(pos.x());
	    gamma_vertex_y.push_back(pos.y());
	    gamma_vertex_z.push_back(pos.z());
	    //std::cout << "OM position: " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;
	    num_om_gamma.push_back(om_number);
	    gamma_type.push_back(it_hit->get_geom_id().get_type());
	  time_gamma.push_back(it_hit->get_time());
	  nb_gamma++;

	  if (event.has("SD")){
	    const snemo::datamodel::unified_digitized_data & UDD = event.get<snemo::datamodel::unified_digitized_data>("UDD");
	    const auto& calorimeter_hits = UDD.get_calorimeter_hits();
	    for (const auto& hit : calorimeter_hits) {
	      if(om_number==snemo::datamodel::om_num(hit->get_geom_id())){
		const datatools::properties& aux = hit->get_auxiliaries();
		for (const auto& key : aux.keys()) {
		  const auto& value = aux.get(key);
		  double real_value = 0.0;
		  if(key=="evis"){
		    value.get_value(real_value);
		    gamma_evis.push_back(real_value);
		  }		 
		  else if(key=="edep"){
		    value.get_value(real_value);
		    gamma_edep.push_back(real_value);
		  }
		  else if(key=="edep_bcu"){
		    value.get_value(real_value);
		    gamma_edep_bcu.push_back(real_value);
		  }		  
		  else if(key=="evis_bc"){
		    value.get_value(real_value);
		    gamma_evis_bc.push_back(real_value);
		  }
		  else if(key=="evis_bcu"){
		    value.get_value(real_value);
		    gamma_evis_bcu.push_back(real_value);
		  }
		  else if(key=="evis_u"){
		    value.get_value(real_value);
		    gamma_evis_u.push_back(real_value);
		  }
		}
		break;
	      }
	    }
	  }
	  }
	}
    
     
      //electron start
      //vector<vector<double>> vertex_gamma_tot;
    
      vector<int> cluster_id;
      for (const datatools::handle<snemo::datamodel::particle_track> & particle : PTD.particles())
	{	 
	  if(particle->get_charge() == snemo::datamodel::particle_track::CHARGE_UNDEFINED || particle->get_charge() == snemo::datamodel::particle_track::CHARGE_POSITIVE || particle->get_charge() == snemo::datamodel::particle_track::CHARGE_NEGATIVE){
	    if(particle->get_charge() == snemo::datamodel::particle_track::CHARGE_POSITIVE){
	      particle_type.push_back(1);
	    }
	    if(particle->get_charge() == snemo::datamodel::particle_track::CHARGE_NEGATIVE){
	      particle_type.push_back(-1);
	    }
	    bool vertex_close_to_the_source = false;
	    bool vertex_associated_to_a_calo = false;
            double x_calo, y_calo, z_calo, x_foil, y_foil, z_foil;	  
	      
	    for(const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices()){
	      if(vertex->is_on_reference_source_plane()){
	      //if(vertex->is_on_calibration_source()){
		x_foil = vertex->get_spot().get_position().getX();
                y_foil = vertex->get_spot().get_position().getY();
                z_foil = vertex->get_spot().get_position().getZ();
		vertex_close_to_the_source=1;
	      }
	      else if(vertex->is_on_source_foil() && vertex_close_to_the_source==0){
		x_foil = vertex->get_spot().get_position().getX();
                y_foil = vertex->get_spot().get_position().getY();
                z_foil = vertex->get_spot().get_position().getZ();
                vertex_close_to_the_source=1;
	      }
	       
	      else if((vertex->is_on_main_calorimeter() || vertex->is_on_x_calorimeter()) && particle->get_associated_calorimeter_hits().size()==1) //we forced MW only analysis
		{
		  if(vertex->is_on_main_calorimeter())
		    {
		      type_elec.push_back("MW");
		    }
		  else
		    {
		      type_elec.push_back("XW");
		    }
		  x_calo = vertex->get_spot().get_position().getX();
		  y_calo = vertex->get_spot().get_position().getY();
		  z_calo = vertex->get_spot().get_position().getZ();		
		  vertex_associated_to_a_calo = true;
		}
	    }
	    int source_num;
	    if(vertex_close_to_the_source==1 && vertex_associated_to_a_calo==1 /*&& particle->get_trajectory_handle()->get_pattern().number_of_kinks()==0*/ /* && compute_ellipse(y_foil,z_foil,source_num)<1 */){
	      double dy, dz;
	      ellipse_source.push_back(compute_ellipse(y_foil,z_foil,source_num,dy,dz));
	      delta_y_source.push_back(dy);
	      delta_z_source.push_back(dz);
	      number_of_kinks_per_track.push_back(particle->get_trajectory_handle()->get_pattern().number_of_kinks());
	      //verify that the cluster from the 2 electrons are not the sames
	      const auto& calorimeter_hits = particle->get_associated_calorimeter_hits();
	      //if(std::find(num_om.begin(), num_om.end(), snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id())) == num_om.end()){
		vertex_3D_start_x.push_back(x_foil);
		vertex_3D_start_y.push_back(y_foil);
		vertex_3D_start_z.push_back(z_foil);
		vertex_3D_end_x.push_back(x_calo);
		vertex_3D_end_y.push_back(y_calo);
		vertex_3D_end_z.push_back(z_calo);
		cluster_id.push_back(particle->get_trajectory().get_cluster().get_hit_id());
		//energy_elec.push_back(calorimeter_hits[0]->get_energy());
		time_elec.push_back(calorimeter_hits[0]->get_time());
		//cout<<calorimeter_hits[0]->get_time()<<endl;
		side_elec.push_back(calorimeter_hits[0]->get_geom_id().get(1));
		int number = snemo::datamodel::om_num(calorimeter_hits[0]->get_geom_id());
	    
		num_om.push_back(number);
		num_source.push_back(source_num);
		delta_y_calo_center.push_back(abs(y_calo - y_map[number]));
		delta_z_calo_center.push_back(abs(z_calo - z_map[number]));
		number_of_electrons++;
		const datatools::handle<snemo::datamodel::precalibrated_calorimeter_hit> & pcd_calo_hit = pCD.calorimeter_hits()[calorimeter_hits[0]->get_hit_id()];
		charge_elec.push_back((-pcd_calo_hit->get_charge())/(1E-9*CLHEP::volt*CLHEP::second * 0.00023841858));
		//charge to nv*s -> 1 unity of UDD charge = 1 unity calo time * 1 unity calo amplitude
		// 0.390625 ns x 0.61035156 mV = 0.00023841858 nVs
		if (event.has("SD")){
		  const snemo::datamodel::unified_digitized_data & UDD = event.get<snemo::datamodel::unified_digitized_data>("UDD");
		  const auto& calorimeter_hits = UDD.get_calorimeter_hits();
		  for (const auto& hit : calorimeter_hits) {
		    if(number==snemo::datamodel::om_num(hit->get_geom_id())){
		      const datatools::properties& aux = hit->get_auxiliaries();
		      for (const auto& key : aux.keys()) {
			const auto& value = aux.get(key);
			if(key=="evis"){
			  double real_value = 0.0;
			  value.get_value(real_value);
			  evis.push_back(real_value);
			}
			else if(key=="evis_bc"){
			  double real_value = 0.0;
			  value.get_value(real_value);
			  evis_bc.push_back(real_value);
			}
			else if(key=="evis_bcu"){
			  double real_value = 0.0;
			  value.get_value(real_value);
			  evis_bcu.push_back(real_value);
			}
			else if(key=="evis_u"){
			  double real_value = 0.0;
			  value.get_value(real_value);
			  evis_u.push_back(real_value);
			}
			else if(key=="edep_bcu"){
			  double real_value = 0.0;
			  value.get_value(real_value);
			  edep_bcu.push_back(real_value);
			}
		      }		      
		      break;		      
		      //}//if number==
		    }
		}
	      }
	    }
	    else{
              cellules_non_associated++;	      
	    }
	  }//end e-	  	  
	  
	  else{
	    cellules_non_associated++;
	  }
	}//end loop on particle
    
  bool evt_save=false;
      
      nb_elec_linked=0;
      for(int i=0; i<number_of_electrons; i++){
	//electron selection from the source
	if(ellipse_source[i]>2) continue;
	nb_elec_linked_vec.push_back(0);
	nb_gamma_per_elec=0;
	if(number_of_kinks_per_track[i]!=0) continue;
	double x_start = vertex_3D_start_x[i];
	double y_start = vertex_3D_start_y[i];
	double z_start = vertex_3D_start_z[i];
	double x_end = vertex_3D_end_x[i];
	double y_end = vertex_3D_end_y[i];
	double z_end = vertex_3D_end_z[i];	
	double e_elec;
	if (event.has("SD")){
	  e_elec = edep_bcu[i];
	}
	else{
	  e_elec = (0.23841858*charge_elec[i] * calib_param[num_om[i]])/1000;
	}
	total_energy_electron.push_back(e_elec);

 	//cout<<charge_elec[i]<<" "<<calib_param[num_om[i]]<<" "<<time_elec[i]<<endl;
	double t_elec = time_elec[i]/ CLHEP::ns;
	double new_c = 299.792458;//change c unit to mm/ns
	double d_e = std::sqrt(pow(x_end-x_start ,2) + pow(y_end - y_start,2) + pow(z_end - z_start,2));
	double e_elec_tot = e_elec+0.511;
	double beta_e = std::sqrt(1.0 - (0.511*0.511)/(e_elec_tot*e_elec_tot));
	//out<<beta_e<<" "<<e_elec_tot<<" "<<e_elec<<endl;
	double v_e = beta_e * new_c;
	double th_time_e = d_e/v_e;
	double th_time_e_c = d_e/new_c;
	double t_abs_e = t_elec-th_time_e;
	double t_abs_e_c = t_elec-th_time_e_c;
	for(int j=0; j<nb_gamma; j++){
	  nb_gamma_per_elec_vec.push_back(0);
	  double e_gamma;
	  if (event.has("SD")){
	    e_gamma = gamma_edep[j];
	    total_energy_gamma_vis.push_back(gamma_evis[j]);	   
	  }
	  else{	    
	    e_gamma = (0.23841858*charge_gamma[j] * calib_param[num_om_gamma[j]])/1000;
	  }
	  num_om_gamma_f.push_back(num_om_gamma[j]);
	  total_energy_gamma.push_back(e_gamma);
	  double x_gamma = gamma_vertex_x[j];
	  double y_gamma = gamma_vertex_y[j];
	  double z_gamma = gamma_vertex_z[j];
	  double t_gamma = time_gamma[j]/ CLHEP::ns;
	  double d_gamma = std::sqrt(pow(x_gamma-x_start,2)+pow(y_gamma-y_start,2)+pow(z_gamma-z_start,2));
	  double t_th_gamma = d_gamma / 299.792458;
	  double t_abs_gamma = t_gamma-t_th_gamma;
	  dt_gamma_elec_c.push_back(t_abs_e_c-t_abs_gamma);
	  base_dt_gamma_elec.push_back(t_elec - t_gamma);
	  if(calib_param[num_om[i]] != 0){
	    dt_gamma_elec.push_back(t_abs_e-t_abs_gamma);
	  }
	  if(!((e_elec>0.7 && e_elec<1.2) || (e_elec>0.3 && e_elec<0.65))) continue;
	  //if(!((e_gamma>0.7 && e_gamma<1.2) || (e_gamma>0.25 && e_gamma<0.65)) || (e_gamma>0.25 && e_gamma<0.65 && e_elec>0.25 && e_elec<0.65)) continue;
	  if(abs(t_abs_e-t_abs_gamma)<2){
	    evt_save=true;
	    nb_gamma_per_elec_vec[j] ++;
	    nb_gamma_per_elec++;
	    energy_gamma.push_back(e_gamma);
	    energy_elec.push_back(e_elec);
	    if (event.has("SD")){
	      gamma_evis_vec.push_back(gamma_evis[j]);
	      gamma_edep_vec.push_back(gamma_edep[j]);
	      gamma_edep_bcu_vec.push_back(gamma_edep_bcu[j]);
	      gamma_evis_bcu_vec.push_back(gamma_evis_bcu[j]);
	      gamma_evis_bc_vec.push_back(gamma_evis_bc[j]);
	      gamma_evis_u_vec.push_back(gamma_evis_u[j]);
	      
	    }
	    else{
	      charge_gamma_finale.push_back(charge_gamma[j]);
	      charge_elec_finale.push_back(charge_elec[i]);
	    }
	  }
	}//for j gamma
	
	if(nb_gamma_per_elec>0){
	  nb_elec_linked_vec[i] ++;
	}
      }//for i elec
      


	//if(number_of_electrons>0){
      if(evt_save==1){
	tree->Fill();
      }
      //}
      //cout<<event_number<<" "<<number_of_electrons<<endl;
      event_number++;
    }
  
  return dpp::base_module::PROCESS_SUCCESS;
}

