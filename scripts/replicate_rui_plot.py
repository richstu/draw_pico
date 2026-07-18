do_vbf = True #False = ggF

#cell 1
if True:
  # when applying scores
  import xgboost as xgb
  from sklearn.metrics import roc_auc_score
  import uproot
  import pandas as pd
  import numpy as np
  import awkward as ak
  import os
  import glob
  import joblib

def concatenate_flat(nparrays: list) -> np.array:
  """Concatenates flat numpy arrays by making new axis
  """
  new_arrays = [arr[:,np.newaxis] for arr in nparrays]
  return np.concatenate(new_arrays, axis=1)

#cell 2
if True:
  # when plotting
  import matplotlib
  import matplotlib.pyplot as plt
  import uproot
  import pandas as pd
  import matplotlib.pylab as pl
  import math
  #from matplotlib import pyplot as plt
  from matplotlib.colors import ListedColormap
  import matplotlib.patches as mpatches
  import matplotlib.markers as mmark
  
  red_patch = mpatches.Patch(color='red', label='Signal')
  green_patch = mpatches.Patch(color='green', label='VBF')
  blue_patch = mpatches.Patch(color='#3f90da', label='SM Z$\\gamma$')
  orange_patch = mpatches.Patch(color='#ffa90e', label='DY')
  teal_patch = mpatches.Patch(color='#92dadd',label='VBS Z')
  from matplotlib.lines import Line2D
  from matplotlib.patches import Rectangle
  import matplotlib.lines as mlines
  #dot = mlines.Line2D([], [], color='black', marker='.', linestyle='None',
  #                          markersize=10, label='Data')
  
  l = Line2D([], [], color='red', lw = 3)
  lg = Line2D([], [], color='green', lw = 3)
  extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
  
  import mplhep as hep
  hep.style.use("CMS")
  plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = False
  #plt.rcParams['tick.left'] = plt.rcParams['ytick.labelleft'] = False
  plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False
  
  import matplotlib.gridspec as gridspec
  from scipy.stats import gamma
  def poisson_interval(n, conf_level=0.6827):
      """
      Reproduce ROOT::kPoisson (Garwood) intervals for counts n.
      n can be scalar or numpy array.
      """
      n = np.asarray(n)
      alpha = 1.0 - conf_level
      
      lower = np.where(
          n > 0,
          gamma.ppf(alpha / 2.0, n),  # shape = n, scale = 1
          0.0
      )
      upper = gamma.isf(alpha / 2.0, n + 1)  # isf = 1 - cdf, same as gamma_quantile_c
      
      return lower, upper
  
  
#cell 3
if False:
  # === SETTINGS ===
  #model_path = "BDT_zero_to_one_jet_0.h5"
  #input_root = "GGF_2016APV_pinnacles_all.root"
  #output_root = "output_with_score.root"
  year_map = [
      ("2016APV", 0),
      ("2016", 1),
      ("2017", 2),
      ("2018", 3),
      ("2022EE", 5),
      ("2022", 4),
      ("2023BPix", 7),
      ("2023", 6),
  ]
  input_folder = "../../"
  #input_folder = "input_all/"
  #output_folder = "Output_ggF_xgb_commonparam/"
  output_folder = "Output_ggF_rui_commonparam_era/"
  
  
  tree_name = "TreeS"  # or "TreeB"
  #feature_order = [
  #"cosTheta", "phi", "costheta",
  #                  "l1_rapidity", "l2_rapidity", "photon_rapidity",
  #                              "min_dR", "max_dR",
  #                               "photon_mva", "photon_res",
  #                               "pt_mass", "llphoton_dijet_dphi", "j1_m",
  #                                "j1_eta", "llphoton_dijet_balance", "njet",
  #                                "photon_zeppenfeld",hadd si
  #           "photon_jet1_dr", "j1_pt", "llphoton_hmiss_photon_dphi"]
  
  
  
  feature_order = [\
    "costheta", \
    "cosTheta", \
    "pt_mass", \
    "l1_rapidity",\
    "l2_rapidity",\
    "photon_rapidity",\
    "phi",\
    "photon_mva",\
    "photon_res",\
    "min_dR",\
    "max_dR",\
    "llphoton_hmiss_photon_dphi",\
    "photon_jet1_dr",\
    "photon_zeppenfeld",\
    "llphoton_dijet_dphi",\
    "llphoton_dijet_balance",\
    "j1_eta",\
    "j1_m",\
    "j1_pt",\
     "njet" , "year"] 
  #h5_models = {
  #    0: "BDT_zero_to_one_jet_0_commonparam.h5",
  #    1: "BDT_zero_to_one_jet_1_commonparam.h5",
  #    2: "BDT_zero_to_one_jet_2_commonparam.h5",
  #    3: "BDT_zero_to_one_jet_3_commonparam.h5",
  #}
  
  h5_models = {
      0: "model_xgb_0_ggf_era.json",
      1: "model_xgb_1_ggf_era.json",
      2: "model_xgb_2_ggf_era.json",
      3: "model_xgb_3_ggf_era.json",
  }
  
  root_files = glob.glob(os.path.join(input_folder, "data_20*_ggf_fixedjet.root"))
  #root_files = glob.glob(os.path.join(input_folder, "*F_20*_all_fixedjet.root"))
  
  for input_root in root_files:
      filename = os.path.basename(input_root)
      
      year_label = None
      for year_str, label in year_map:
          if year_str in filename:
              year_label = label
              break
  
      if year_label is None:
          print(f"[WARNING] No year label found for {filename}, skipping.")
          continue
      
      output_filename = filename.replace("_pinnacles_ggf_fixedjet", "_output")
  #    output_filename = filename.replace("_pinnacles_all_fixedjet", "_output")
      output_root = os.path.join(output_folder, output_filename)
      print(f"Processing {filename} → {output_root} (year={year_label})")
      
  #    continue
      # Detect tree and output folder
  # === 2. Load ROOT file and tree ===
      with uproot.open(input_root) as file:
          if "TreeS" in file:
              tree = file["TreeS"]
          elif "TreeB" in file:
              tree = file["TreeB"]
      # Load all branches for output, but only feature_order for prediction
          full_df = tree.arrays(library="pd")
          # Add year column
          full_df["year"] = year_label
          
  #    feature_df = full_df[feature_order]
  #    full_df = tree.arrays(library="pd")  # full_df contains all branches
  
  # === 2. Safety check ===
      assert 'part' in full_df.columns, "'part' branch is missing!"
  
  # === 3. Create BDT_score column ===
      full_df["BDT_score"] = np.nan  # Initialize with NaNs   
      full_df["BDT_score_val"] = np.nan    
      full_df["BDT_score_train"] = np.nan    
      full_df["BDT_score_train_1"] = np.nan 
      
   #   full_df["BDT_score_t"] = np.nan  # Initialize with NaNs   
   #   full_df["BDT_score_val_t"] = np.nan    
   #   full_df["BDT_score_train_t"] = np.nan    
   #   full_df["BDT_score_train_1_t"] = np.nan 
      
  
      for fold, model_path in h5_models.items():
          print(f"Applying model for fold {fold}: {model_path}")
      
      # Load model
          booster = xgb.Booster()
          booster.load_model(model_path)
      
      # Select subset for this fold
  #    feature_df = full_df.loc[fold_mask, feature_order]  # just features
  #        transformer = joblib.load(f"BDT_tsf_zero_to_one_jet_{fold}_commonparam.pkl")
          part = full_df["part"]
          cut = (full_df["njet"] < 2) & (full_df["met"] < 90)
          test_mask = (part == fold) & cut
          val_mask = (part == (fold + 1) % 4) & cut
          train_mask =  (part == (fold + 2) % 4) & cut
          train_mask_1 = (part == (fold + 3) % 4) & cut
          
  #        column1= ["BDT_score_t","BDT_score_val_t","BDT_score_train_t","BDT_score_train_1_t"]
          i = 0
      # Apply only to events for this fold
          for mask, column in [(test_mask, "BDT_score"),
                               (val_mask, "BDT_score_val"),
                               (train_mask, "BDT_score_train"), 
                               (train_mask_1, "BDT_score_train_1")]:
              if mask.any():
                  dmatrix = xgb.DMatrix(full_df.loc[mask, feature_order])
                  preds = booster.predict(dmatrix)
                  full_df.loc[mask, column] = preds
    #              preds_transformed = transformer.transform(preds.reshape(-1, 1))
    #              full_df.loc[mask, column1[i]] = preds_transformed.ravel()
                  i += 1
                  #xgb_scores = model.predict(dmatrix)  # shape (n_samples,)
  
  
  # === 4. Save to ROOT
      ak_array = ak.Array(full_df.to_dict(orient="list"))
      with uproot.recreate(output_root) as fout:
          fout["outtree"] = ak_array
  
      
  print(f"✅ Done! BDT_score written to {output_root} in tree {tree_name}")

#cell 4
if False:
  #input_file = "input_all/SM1_2022EE_pinnacles_all.root"
  input_file = "../../SM1_2022_pinnacles_ggf_fixedjet.root"
  output_file = "Output_ggF_rui_commonparam_era/SM1_2022_output.root"
  
  
  year_label = None
  for year_str, label in year_map:
      if year_str in input_file:
          year_label = label
          break
  
  if year_label is None:
      print(f"[WARNING] No year label found for {input_file}, skipping.")
  
  
  #input_file  = "../../DYmix_2018_pinnacles_ggf_fixedjet_v1.root"
  #output_file = "Output_ggF_rui_commonparam/DYmix_2018_output.root" 
  tree_name = "TreeB"  # or TreeB
  #feature_order = [
  #"cosTheta", "phi", "costheta",
  #                  "l1_rapidity", "l2_rapidity", "photon_rapidity",
  #                              "min_dR", "max_dR",
  #                               "photon_mvann", "photon_res",
  #                               "pt_mass", "llphoton_dijet_dphi", "j1_m",
  #                                "j1_eta", "llphoton_dijet_balance", "njet",
  #                                "dijet_dphi", "photon_zeppenfeld",
  #           "photon_jet1_dr", "j1_pt", "llphoton_hmiss_photon_dphi"]
  #feature_order = [
  #"cosTheta", "phi", "costheta",
  #                  "l1_rapidity", "l2_rapidity", "photon_rapidity",
  #                              "min_dR", "max_dR",
  #                               "photon_mva", "photon_res",
  #                               "pt_mass", "llphoton_dijet_dphi", "j1_m",
  #                                "j1_eta", "llphoton_dijet_balance", "njet",
  #                                "photon_zeppenfeld",
  #           "photon_jet1_dr", "j1_pt", "llphoton_hmiss_photon_dphi"]
  
  feature_order = [\
      "costheta", \
    "cosTheta", \
    "pt_mass", \
    "l1_rapidity",\
    "l2_rapidity",\
    "photon_rapidity",\
    "phi",\
    "photon_mva",\
    "photon_res",\
    "min_dR",\
    "max_dR",\
    "llphoton_hmiss_photon_dphi",\
    "photon_jet1_dr",\
    "photon_zeppenfeld",\
    "llphoton_dijet_dphi",\
    "llphoton_dijet_balance",\
    "j1_eta",\
    "j1_m",\
    "j1_pt",\
     "njet" , "year"] 
  #h5_models = {
  #    0: "BDT_zero_to_one_jet_0_commonparam.h5",
  #    1: "BDT_zero_to_one_jet_1_commonparam.h5",
  #    2: "BDT_zero_to_one_jet_2_commonparam.h5",
  #    3: "BDT_zero_to_one_jet_3_commonparam.h5",
  #}
  
  h5_models = {
      0: "model_xgb_0_ggf_era.json",
      1: "model_xgb_1_ggf_era.json",
      2: "model_xgb_2_ggf_era.json",
      3: "model_xgb_3_ggf_era.json",
  }
  
  with uproot.recreate(output_file) as fout:
      tree_initialized = False  # fclag to track if tree is created
  
  
      # === Chunk iterator: load all branches in steps
      for arrays in uproot.iterate(
          {input_file: tree_name},
          expressions=None,  # load all branches
          step_size="100 MB",  # or int number of entries
          library="pd"
      ):
          # Backup full set for output
          full_df = arrays.copy()
          full_df["year"] = year_label
  
          if "part" not in full_df.columns:
              raise RuntimeError("Missing 'part' column in input tree.")
  
  # === 3. Create BDT_score column ===
          full_df["BDT_score"] = np.nan  # Initialize with NaNs   
          full_df["BDT_score_val"] = np.nan    
          full_df["BDT_score_train"] = np.nan    
          full_df["BDT_score_train_1"] = np.nan    
  #        full_df["BDT_score_t"] = np.nan  # Initialize with NaNs   
  #        full_df["BDT_score_val_t"] = np.nan    
  #        full_df["BDT_score_train_t"] = np.nan    
  #        full_df["BDT_score_train_1_t"] = np.nan 
  
          for fold, model_path in h5_models.items():
              print(f"Applying model for fold {fold}: {model_path}")
      
      # Load model
              booster = xgb.Booster()
              booster.load_model(model_path)
      
      # Select subset for this fold
  #    feature_df = full_df.loc[fold_mask, feature_order]  # just features
  #            transformer = joblib.load(f"BDT_tsf_zero_to_one_jet_{fold}_commonparam.pkl")
              part = full_df["part"]
              cut = (full_df["njet"] < 2) & (full_df["met"] < 90)
              test_mask = (part == fold) & cut
              val_mask = (part == (fold + 1) % 4) & cut
              train_mask =  (part == (fold + 2) % 4) & cut
              train_mask_1 = (part == (fold + 3) % 4) & cut
  #            column1= ["BDT_score_t","BDT_score_val_t","BDT_score_train_t","BDT_score_train_1_t"]
              i = 0
  
      # Apply only to events for this fold
              for mask, column in [(test_mask, "BDT_score"),
                                   (val_mask, "BDT_score_val"),
                                   (train_mask, "BDT_score_train"), 
                                   (train_mask_1, "BDT_score_train_1")]:
                  if mask.any():
                      dmatrix = xgb.DMatrix(full_df.loc[mask, feature_order])
                      preds = booster.predict(dmatrix)
                      full_df.loc[mask, column] = preds
  #                    preds_transformed = transformer.transform(preds.reshape(-1, 1))
  #                    full_df.loc[mask, column1[i]] = preds_transformed.ravel()
                      i += 1
   # Write this chunk to the output file
          out_array = ak.Array(full_df.to_dict(orient="list"))
          if not tree_initialized:
              fout.mktree("outtree", {key: out_array[key].type for key in out_array.fields})
              tree_initialized = True
          fout["outtree"].extend(out_array)
      
  
  print(f"✅ Done! BDT_score written to {output_file} in tree outtree")

#cell 5
if False:
    "#input_file = \"input_all/SM1_2022EE_pinnacles_all.root\"\n",
    "input_file = \"../../SM1_2022_pinnacles_ggf_fixedjet.root\"\n",
    "output_file = \"Output_ggF_rui_commonparam_era/SM1_2022_output.root\"\n",
    "\n",
    "\n",
    "year_label = None\n",
    "for year_str, label in year_map:\n",
    "    if year_str in input_file:\n",
    "        year_label = label\n",
    "        break\n",
    "\n",
    "if year_label is None:\n",
    "    print(f\"[WARNING] No year label found for {input_file}, skipping.\")\n",
    "\n",
    "\n",
    "#input_file  = \"../../DYmix_2018_pinnacles_ggf_fixedjet_v1.root\"\n",
    "#output_file = \"Output_ggF_rui_commonparam/DYmix_2018_output.root\" \n",
    "tree_name = \"TreeB\"  # or TreeB\n",
    "#feature_order = [\n",
    "#\"cosTheta\", \"phi\", \"costheta\",\n",
    "#                  \"l1_rapidity\", \"l2_rapidity\", \"photon_rapidity\",\n",
    "#                              \"min_dR\", \"max_dR\",\n",
    "#                               \"photon_mvann\", \"photon_res\",\n",
    "#                               \"pt_mass\", \"llphoton_dijet_dphi\", \"j1_m\",\n",
    "#                                \"j1_eta\", \"llphoton_dijet_balance\", \"njet\",\n",
    "#                                \"dijet_dphi\", \"photon_zeppenfeld\",\n",
    "#           \"photon_jet1_dr\", \"j1_pt\", \"llphoton_hmiss_photon_dphi\"]\n",
    "#feature_order = [\n",
    "#\"cosTheta\", \"phi\", \"costheta\",\n",
    "#                  \"l1_rapidity\", \"l2_rapidity\", \"photon_rapidity\",\n",
    "#                              \"min_dR\", \"max_dR\",\n",
    "#                               \"photon_mva\", \"photon_res\",\n",
    "#                               \"pt_mass\", \"llphoton_dijet_dphi\", \"j1_m\",\n",
    "#                                \"j1_eta\", \"llphoton_dijet_balance\", \"njet\",\n",
    "#                                \"photon_zeppenfeld\",\n",
    "#           \"photon_jet1_dr\", \"j1_pt\", \"llphoton_hmiss_photon_dphi\"]\n",
    "\n",
    "feature_order = [\\\n",
    "    \"costheta\", \\\n",
    "  \"cosTheta\", \\\n",
    "  \"pt_mass\", \\\n",
    "  \"l1_rapidity\",\\\n",
    "  \"l2_rapidity\",\\\n",
    "  \"photon_rapidity\",\\\n",
    "  \"phi\",\\\n",
    "  \"photon_mva\",\\\n",
    "  \"photon_res\",\\\n",
    "  \"min_dR\",\\\n",
    "  \"max_dR\",\\\n",
    "  \"llphoton_hmiss_photon_dphi\",\\\n",
    "  \"photon_jet1_dr\",\\\n",
    "  \"photon_zeppenfeld\",\\\n",
    "  \"llphoton_dijet_dphi\",\\\n",
    "  \"llphoton_dijet_balance\",\\\n",
    "  \"j1_eta\",\\\n",
    "  \"j1_m\",\\\n",
    "  \"j1_pt\",\\\n",
    "   \"njet\" , \"year\"] \n",
    "#h5_models = {\n",
    "#    0: \"BDT_zero_to_one_jet_0_commonparam.h5\",\n",
    "#    1: \"BDT_zero_to_one_jet_1_commonparam.h5\",\n",
    "#    2: \"BDT_zero_to_one_jet_2_commonparam.h5\",\n",
    "#    3: \"BDT_zero_to_one_jet_3_commonparam.h5\",\n",
    "#}\n",
    "\n",
    "h5_models = {\n",
    "    0: \"model_xgb_0_ggf_era.json\",\n",
    "    1: \"model_xgb_1_ggf_era.json\",\n",
    "    2: \"model_xgb_2_ggf_era.json\",\n",
    "    3: \"model_xgb_3_ggf_era.json\",\n",
    "}\n",
    "\n",
    "with uproot.recreate(output_file) as fout:\n",
    "    tree_initialized = False  # fclag to track if tree is created\n",
    "\n",
    "\n",
    "    # === Chunk iterator: load all branches in steps\n",
    "    for arrays in uproot.iterate(\n",
    "        {input_file: tree_name},\n",
    "        expressions=None,  # load all branches\n",
    "        step_size=\"100 MB\",  # or int number of entries\n",
    "        library=\"pd\"\n",
    "    ):\n",
    "        # Backup full set for output\n",
    "        full_df = arrays.copy()\n",
    "        full_df[\"year\"] = year_label\n",
    "\n",
    "        if \"part\" not in full_df.columns:\n",
    "            raise RuntimeError(\"Missing 'part' column in input tree.\")\n",
    "\n",
    "# === 3. Create BDT_score column ===\n",
    "        full_df[\"BDT_score\"] = np.nan  # Initialize with NaNs   \n",
    "        full_df[\"BDT_score_val\"] = np.nan    \n",
    "        full_df[\"BDT_score_train\"] = np.nan    \n",
    "        full_df[\"BDT_score_train_1\"] = np.nan    \n",
    "#        full_df[\"BDT_score_t\"] = np.nan  # Initialize with NaNs   \n",
    "#        full_df[\"BDT_score_val_t\"] = np.nan    \n",
    "#        full_df[\"BDT_score_train_t\"] = np.nan    \n",
    "#        full_df[\"BDT_score_train_1_t\"] = np.nan \n",
    "\n",
    "        for fold, model_path in h5_models.items():\n",
    "            print(f\"Applying model for fold {fold}: {model_path}\")\n",
    "    \n",
    "    # Load model\n",
    "            booster = xgb.Booster()\n",
    "            booster.load_model(model_path)\n",
    "    \n",
    "    # Select subset for this fold\n",
    "#    feature_df = full_df.loc[fold_mask, feature_order]  # just features\n",
    "#            transformer = joblib.load(f\"BDT_tsf_zero_to_one_jet_{fold}_commonparam.pkl\")\n",
    "            part = full_df[\"part\"]\n",
    "            cut = (full_df[\"njet\"] < 2) & (full_df[\"met\"] < 90)\n",
    "            test_mask = (part == fold) & cut\n",
    "            val_mask = (part == (fold + 1) % 4) & cut\n",
    "            train_mask =  (part == (fold + 2) % 4) & cut\n",
    "            train_mask_1 = (part == (fold + 3) % 4) & cut\n",
    "#            column1= [\"BDT_score_t\",\"BDT_score_val_t\",\"BDT_score_train_t\",\"BDT_score_train_1_t\"]\n",
    "            i = 0\n",
    "\n",
    "    # Apply only to events for this fold\n",
    "            for mask, column in [(test_mask, \"BDT_score\"),\n",
    "                                 (val_mask, \"BDT_score_val\"),\n",
    "                                 (train_mask, \"BDT_score_train\"), \n",
    "                                 (train_mask_1, \"BDT_score_train_1\")]:\n",
    "                if mask.any():\n",
    "                    dmatrix = xgb.DMatrix(full_df.loc[mask, feature_order])\n",
    "                    preds = booster.predict(dmatrix)\n",
    "                    full_df.loc[mask, column] = preds\n",
    "#                    preds_transformed = transformer.transform(preds.reshape(-1, 1))\n",
    "#                    full_df.loc[mask, column1[i]] = preds_transformed.ravel()\n",
    "                    i += 1\n",
    " # Write this chunk to the output file\n",
    "        out_array = ak.Array(full_df.to_dict(orient=\"list\"))\n",
    "        if not tree_initialized:\n",
    "            fout.mktree(\"outtree\", {key: out_array[key].type for key in out_array.fields})\n",
    "            tree_initialized = True\n",
    "        fout[\"outtree\"].extend(out_array)\n",
    "    \n",
    "\n",
    "print(f\"✅ Done! BDT_score written to {output_file} in tree outtree\")\n"

#cell ?
if True:
    #dir_ = "Output_ggF_xgb_fixedjet_neweight/"
    
    #dir_ = "Output_ggF_rui_commonparam/"
    #dir_1 = "Output_ggF_rui_redwood_v1_val/"
    #dir_ = "Output_ggF_rui_redwood_v2/"
    #=========================================================================
    # MO modification
    #=========================================================================
    dir_1 = "~rz393/workspace_619/h/Output_ggF_rui_redwood_v1_val/"
    dir_ = "~rz393/workspace_619/h/Output_ggF_rui_redwood_v2/"
    if do_vbf:
      dir_1 = "Output_VBF_rui_redwood_v1_ext_val/"
      dir_ = "~rz393/workspace_619/h/Output_VBF_rui_redwood_v2/"
    #=========================================================================
    #dir_ = "Output_ggF_rui_commonparam_redwood"
    #dir_ = "Output_ggF_rui_commonparam_corr_r2_redwood/"
    #dir_ = "Output_VBF_xgb_neweight/"
    #dir_xgb = "relpt_peking_new"
    #dir_xgb = "relpt_peking_run2p3"
    
    #dir_ = "../../xgboost_pku/"
    #dir_xgb = "output_relpt"
      
    variables = ["BDT_score_val","BDT_score","BDT_score_train","BDT_score_train_1",#"BDT_score_val_t","BDT_score_t","BDT_score_train_t","BDT_score_train_1_t",\
                 "llphoton_m","llphoton_refit_m",\
                 "weight_corr","met", "part", "ll_lepid", "njet", "nbdfm", "photon_jet1_dr", "nel","nmu"]#, "weight"]
    
    
    
    arr_GGF_xgb_ggf = uproot.concatenate(dir_+"/GGF_20*_output.root:outtree", variables)
    arr_GGF_xgb_ggf["index"] = 1
    arr_GGF_xgb_ggf["weight1"] = arr_GGF_xgb_ggf["weight_corr"]#*arr_GGF_xgb["w_year"]
    arr_GGF_xgb_ggf["type"] = 0
    arr_GGF_xgb_ggf["lep"] = arr_GGF_xgb_ggf["ll_lepid"]
    
    arr_VBF_xgb_ggf = uproot.concatenate(dir_+"/VBF_20*_output.root:outtree", variables)
    arr_VBF_xgb_ggf["weight1"] = arr_VBF_xgb_ggf["weight_corr"]#*arr_VBF_xgb["w_year"]
    arr_VBF_xgb_ggf["index"] = 1
    arr_VBF_xgb_ggf["type"] = 1
    arr_VBF_xgb_ggf["lep"] = arr_VBF_xgb_ggf["ll_lepid"]
    
    
    arr_SMZg_xgb_ggf = uproot.concatenate(dir_+"/SM*_20*_output.root:outtree", variables)
    arr_SMZg_xgb_ggf["weight1"] = arr_SMZg_xgb_ggf["weight_corr"]#*arr_SMZg_xgb["w_year"]
    arr_SMZg_xgb_ggf["index"] = 0
    arr_SMZg_xgb_ggf["type"] = 0
    arr_SMZg_xgb_ggf["lep"] = arr_SMZg_xgb_ggf["ll_lepid"]
    
    #arr_DY_xgb_ggf = uproot.concatenate(dir_+"/DY0_20*_output.root:outtree", variables)
    suffix = '_1'
    #if do_vbf:
    #    suffix = ''
    arr_DY_xgb_ggf = uproot.concatenate(dir_1+"/DY*_20*_output"+suffix+".root:outtree", variables)
    arr_DY_xgb_ggf["weight1"] = arr_DY_xgb_ggf["weight_corr"]#*arr_DY_xgb["w_year"]
    arr_DY_xgb_ggf["index"] = 0
    arr_DY_xgb_ggf["type"] = 1
    arr_DY_xgb_ggf["lep"] = arr_DY_xgb_ggf["ll_lepid"]
    
    
    arr_EWK_xgb_ggf = uproot.concatenate(dir_+"/EWK_20*_output.root:outtree", variables)
    arr_EWK_xgb_ggf["weight1"] = arr_EWK_xgb_ggf["weight_corr"]#*arr_EWK_xgb["w_year"]
    arr_EWK_xgb_ggf["index"] = 0
    arr_EWK_xgb_ggf["type"] = 2
    arr_EWK_xgb_ggf["lep"] = arr_EWK_xgb_ggf["ll_lepid"]
    
    
    arr_data_xgb_ggf = uproot.concatenate(dir_+"/data_20*_output.root:outtree", variables)
    arr_data_xgb_ggf["weight1"] = 1
    arr_data_xgb_ggf["index"] = 3
    arr_data_xgb_ggf["type"] = 0
    arr_data_xgb_ggf["lep"] = arr_data_xgb_ggf["ll_lepid"]
    #print(arr_VBF["BDT_score_val_t"]["BDT_score_val_t"])
    #arr_data=uproot.concatenate(dir_+"/data_20*_output.root:outtree", output_vars)
    #arr_data["weight1"] = 1
    
    
    mllg_low = 120
    mllg_high = 130
    
     
    mask_GGF_xgb_ggf = (arr_GGF_xgb_ggf["met"] <90) & (arr_GGF_xgb_ggf["njet"]<2) & ((arr_GGF_xgb_ggf["nel"]+arr_GGF_xgb_ggf["nmu"])==2) #& ((arr_GGF_xgb_ggf["njet"] == 0) |(arr_GGF_xgb_ggf["photon_jet1_dr"]> 0.04))
    mask_VBF_xgb_ggf = (arr_VBF_xgb_ggf["met"] <90) & (arr_VBF_xgb_ggf["njet"]<2) & ((arr_VBF_xgb_ggf["nel"]+arr_VBF_xgb_ggf["nmu"])==2)#&  ((arr_VBF_xgb_ggf["njet"] == 0) |(arr_VBF_xgb_ggf["photon_jet1_dr"]> 0.04))
    mask_SMZg_xgb_ggf = (arr_SMZg_xgb_ggf["met"] <90) & (arr_SMZg_xgb_ggf["njet"]<2)& ((arr_SMZg_xgb_ggf["nel"]+arr_SMZg_xgb_ggf["nmu"])==2)
    mask_DY_xgb_ggf = (arr_DY_xgb_ggf["met"] <90) & (arr_DY_xgb_ggf["njet"]<2)& ((arr_DY_xgb_ggf["nel"]+arr_DY_xgb_ggf["nmu"])==2)
    mask_EWK_xgb_ggf = (arr_EWK_xgb_ggf["met"] <90) & (arr_EWK_xgb_ggf["njet"]<2)& ((arr_EWK_xgb_ggf["nel"]+arr_EWK_xgb_ggf["nmu"])==2)
    #mask_data_xgb_ggf = ((arr_data_xgb_ggf["llphoton_refit_m"] <= mllg_low) | (arr_data_xgb_ggf["llphoton_refit_m"] >= mllg_high)) & (arr_data_xgb_ggf["met"] <90) & (arr_data_xgb_ggf["njet"] < 2) & ((arr_data_xgb_ggf["nel"]+arr_data_xgb_ggf["nmu"])==2)#\
    mask_data_xgb_ggf = (arr_data_xgb_ggf["met"] <90) & (arr_data_xgb_ggf["njet"] < 2) & ((arr_data_xgb_ggf["nel"]+arr_data_xgb_ggf["nmu"])==2)#\
    if do_vbf:
      mask_GGF_xgb_ggf = (arr_GGF_xgb_ggf['nbdfm']<1) & (arr_GGF_xgb_ggf["njet"]>=2) & ((arr_GGF_xgb_ggf["nel"]+arr_GGF_xgb_ggf["nmu"])==2)
      mask_VBF_xgb_ggf = (arr_VBF_xgb_ggf['nbdfm']<1) & (arr_VBF_xgb_ggf["njet"]>=2) & ((arr_VBF_xgb_ggf["nel"]+arr_VBF_xgb_ggf["nmu"])==2)
      mask_SMZg_xgb_ggf = (arr_SMZg_xgb_ggf['nbdfm']<1) & (arr_SMZg_xgb_ggf["njet"]>=2)& ((arr_SMZg_xgb_ggf["nel"]+arr_SMZg_xgb_ggf["nmu"])==2)
      mask_DY_xgb_ggf = (arr_DY_xgb_ggf['nbdfm']<1) & (arr_DY_xgb_ggf["njet"]>=2)& ((arr_DY_xgb_ggf["nel"]+arr_DY_xgb_ggf["nmu"])==2)
      mask_EWK_xgb_ggf = (arr_EWK_xgb_ggf['nbdfm']<1) & (arr_EWK_xgb_ggf["njet"]>=2)& ((arr_EWK_xgb_ggf["nel"]+arr_EWK_xgb_ggf["nmu"])==2)
      mask_data_xgb_ggf = (arr_data_xgb_ggf['nbdfm']<1) & (arr_data_xgb_ggf["njet"]>=2) & ((arr_data_xgb_ggf["nel"]+arr_data_xgb_ggf["nmu"])==2)#\
    
    #& ((arr_data_xgb_ggf["njet"] == 0) |(arr_data_xgb_ggf["photon_jet1_dr"]> 0.04))
    
    #mask_GGF_xgb_ggf = (arr_GGF_xgb_ggf["njet"]>=2) &(arr_GGF_xgb_ggf["nbdfm"]==0)
    #mask_VBF_xgb_ggf = (arr_VBF_xgb_ggf["njet"]>=2) &(arr_VBF_xgb_ggf["nbdfm"]==0) 
    #mask_SMZg_xgb_ggf = (arr_SMZg_xgb_ggf["njet"]>=2) &(arr_SMZg_xgb_ggf["nbdfm"]==0)
    #mask_DY_xgb_ggf = (arr_DY_xgb_ggf["njet"]>=2) &(arr_DY_xgb_ggf["nbdfm"]==0)
    #mask_EWK_xgb_ggf = (arr_EWK_xgb_ggf["njet"]>=2) &(arr_EWK_xgb_ggf["nbdfm"]==0)
    #mask_data_xgb_ggf = ((arr_data_xgb_ggf["llphoton_refit_m"] < mllg_low) | (arr_data_xgb_ggf["llphoton_refit_m"] > mllg_high))  & (arr_data_xgb_ggf["njet"] >= 2) & (arr_data_xgb_ggf["nbdfm"] ==0)
    
    
    #print(arr_data_xgb["llphoton_m"][mask_data_xgb])
    
    val = "_val_t"#"_val"
    bdt_score = "BDT_score_"+val
    bdt = "BDT_score_t"
    
        #"twoj":arr_VBF["BDT_score_ggf"+val][mask_VBF],
    
    df_VBF_xgb_ggf= pd.DataFrame({
        "ggf_xgb":arr_VBF_xgb_ggf["BDT_score_val"][mask_VBF_xgb_ggf],
        "ggf_xgb_test":arr_VBF_xgb_ggf["BDT_score"][mask_VBF_xgb_ggf],
    #    "ggf_xgb_t":arr_VBF_xgb_ggf["BDT_score"+val][mask_VBF_xgb_ggf],
    #    "ggf_xgb_test_t":arr_VBF_xgb_ggf["BDT_score_t"][mask_VBF_xgb_ggf],
        "ggf_xgb_train":arr_VBF_xgb_ggf["BDT_score_train"][mask_VBF_xgb_ggf],
        "ggf_xgb_train_1":arr_VBF_xgb_ggf["BDT_score_train_1"][mask_VBF_xgb_ggf],
        "mllg":arr_VBF_xgb_ggf["llphoton_m"][mask_VBF_xgb_ggf],
        "mllg_r":arr_VBF_xgb_ggf["llphoton_refit_m"][mask_VBF_xgb_ggf],
        "weight":arr_VBF_xgb_ggf["weight1"][mask_VBF_xgb_ggf],
        "part":arr_VBF_xgb_ggf["part"][mask_VBF_xgb_ggf],
        "type":arr_VBF_xgb_ggf["type"][mask_VBF_xgb_ggf],
        "lep":arr_VBF_xgb_ggf["lep"][mask_VBF_xgb_ggf]
        }, index=arr_VBF_xgb_ggf["index"][mask_VBF_xgb_ggf])
    
        #"twoj":arr_GGF["BDT_score_ggf"+val][mask_GGF],
    
    df_GGF_xgb_ggf= pd.DataFrame({
        "ggf_xgb":arr_GGF_xgb_ggf["BDT_score_val"][mask_GGF_xgb_ggf],
        "ggf_xgb_test":arr_GGF_xgb_ggf["BDT_score"][mask_GGF_xgb_ggf],
    #    "ggf_xgb_t":arr_GGF_xgb_ggf["BDT_score"+val][mask_GGF_xgb_ggf],
    #    "ggf_xgb_test_t":arr_GGF_xgb_ggf["BDT_score_t"][mask_GGF_xgb_ggf],
        "ggf_xgb_train":arr_GGF_xgb_ggf["BDT_score_train"][mask_GGF_xgb_ggf],
        "ggf_xgb_train_1":arr_GGF_xgb_ggf["BDT_score_train_1"][mask_GGF_xgb_ggf],    
        "mllg":arr_GGF_xgb_ggf["llphoton_m"][mask_GGF_xgb_ggf],
        "mllg_r":arr_GGF_xgb_ggf["llphoton_refit_m"][mask_GGF_xgb_ggf],
        "type":arr_GGF_xgb_ggf["type"][mask_GGF_xgb_ggf],
        "part":arr_GGF_xgb_ggf["part"][mask_GGF_xgb_ggf],
        "lep":arr_GGF_xgb_ggf["lep"][mask_GGF_xgb_ggf],
        "weight":arr_GGF_xgb_ggf["weight1"][mask_GGF_xgb_ggf]
        }, index=arr_GGF_xgb_ggf["index"][mask_GGF_xgb_ggf])
    
    
    df_SMZg_xgb_ggf= pd.DataFrame({
        "ggf_xgb":arr_SMZg_xgb_ggf["BDT_score_val"][mask_SMZg_xgb_ggf],
        "ggf_xgb_test":arr_SMZg_xgb_ggf["BDT_score"][mask_SMZg_xgb_ggf], 
    #    "ggf_xgb_t":arr_SMZg_xgb_ggf["BDT_score"+val][mask_SMZg_xgb_ggf],
    #    "ggf_xgb_test_t":arr_SMZg_xgb_ggf["BDT_score_t"][mask_SMZg_xgb_ggf],
        "ggf_xgb_train":arr_SMZg_xgb_ggf["BDT_score_train"][mask_SMZg_xgb_ggf],
        "ggf_xgb_train_1":arr_SMZg_xgb_ggf["BDT_score_train_1"][mask_SMZg_xgb_ggf], 
        "mllg":arr_SMZg_xgb_ggf["llphoton_m"][mask_SMZg_xgb_ggf],
        "mllg_r":arr_SMZg_xgb_ggf["llphoton_refit_m"][mask_SMZg_xgb_ggf],
        "type":arr_SMZg_xgb_ggf["type"][mask_SMZg_xgb_ggf],
        "part":arr_SMZg_xgb_ggf["part"][mask_SMZg_xgb_ggf],
        "lep":arr_SMZg_xgb_ggf["lep"][mask_SMZg_xgb_ggf],
        "weight":arr_SMZg_xgb_ggf["weight1"][mask_SMZg_xgb_ggf]
        }, index=arr_SMZg_xgb_ggf["index"][mask_SMZg_xgb_ggf])
    
        #"twoj":arr_SMZg["BDT_score_ggf"+val][mask_SMZg],
        
    
    df_DY_xgb_ggf= pd.DataFrame({
        "ggf_xgb":arr_DY_xgb_ggf["BDT_score_val"][mask_DY_xgb_ggf],
        "ggf_xgb_test":arr_DY_xgb_ggf["BDT_score"][mask_DY_xgb_ggf],
    #    "ggf_xgb_t":arr_DY_xgb_ggf["BDT_score"+val][mask_DY_xgb_ggf],
    #    "ggf_xgb_test_t":arr_DY_xgb_ggf["BDT_score_t"][mask_DY_xgb_ggf],
        "ggf_xgb_train":arr_DY_xgb_ggf["BDT_score_train"][mask_DY_xgb_ggf],
        "ggf_xgb_train_1":arr_DY_xgb_ggf["BDT_score_train_1"][mask_DY_xgb_ggf], 
        "mllg":arr_DY_xgb_ggf["llphoton_m"][mask_DY_xgb_ggf],
        "mllg_r":arr_DY_xgb_ggf["llphoton_refit_m"][mask_DY_xgb_ggf],
        "type":arr_DY_xgb_ggf["type"][mask_DY_xgb_ggf],
        "part":arr_DY_xgb_ggf["part"][mask_DY_xgb_ggf],
        "lep":arr_DY_xgb_ggf["lep"][mask_DY_xgb_ggf],
        "weight":arr_DY_xgb_ggf["weight1"][mask_DY_xgb_ggf]
        }, index=arr_DY_xgb_ggf["index"][mask_DY_xgb_ggf])
    
        #"twoj":arr_DY["BDT_score_ggf"+val][mask_DY],
        
    df_EWK_xgb_ggf= pd.DataFrame({
        "ggf_xgb":arr_EWK_xgb_ggf["BDT_score_val"][mask_EWK_xgb_ggf],
        "ggf_xgb_test":arr_EWK_xgb_ggf["BDT_score"][mask_EWK_xgb_ggf], 
    #    "ggf_xgb_t":arr_EWK_xgb_ggf["BDT_score"+val][mask_EWK_xgb_ggf],
    #    "ggf_xgb_test_t":arr_EWK_xgb_ggf["BDT_score_t"][mask_EWK_xgb_ggf], 
        "ggf_xgb_train":arr_EWK_xgb_ggf["BDT_score_train"][mask_EWK_xgb_ggf],
        "ggf_xgb_train_1":arr_EWK_xgb_ggf["BDT_score_train_1"][mask_EWK_xgb_ggf], 
        "mllg":arr_EWK_xgb_ggf["llphoton_m"][mask_EWK_xgb_ggf],
        "mllg_r":arr_EWK_xgb_ggf["llphoton_refit_m"][mask_EWK_xgb_ggf],
        "type":arr_EWK_xgb_ggf["type"][mask_EWK_xgb_ggf],
        "part":arr_EWK_xgb_ggf["part"][mask_EWK_xgb_ggf],
        "lep":arr_EWK_xgb_ggf["lep"][mask_EWK_xgb_ggf],
        "weight":arr_EWK_xgb_ggf["weight1"][mask_EWK_xgb_ggf]
        }, index=arr_EWK_xgb_ggf["index"][mask_EWK_xgb_ggf])
    
    df_data_xgb_ggf= pd.DataFrame({
        "ggf_xgb":arr_data_xgb_ggf["BDT_score_val"][mask_data_xgb_ggf],
        "ggf_xgb_test":arr_data_xgb_ggf["BDT_score"][mask_data_xgb_ggf],
    #    "ggf_xgb_t":arr_data_xgb_ggf["BDT_score"+val][mask_data_xgb_ggf],
    #    "ggf_xgb_test_t":arr_data_xgb_ggf["BDT_score_t"][mask_data_xgb_ggf],
        "ggf_xgb_train":arr_data_xgb_ggf["BDT_score_train"][mask_data_xgb_ggf],
        "ggf_xgb_train_1":arr_data_xgb_ggf["BDT_score_train_1"][mask_data_xgb_ggf],
        "mllg":arr_data_xgb_ggf["llphoton_m"][mask_data_xgb_ggf],
        "type":arr_data_xgb_ggf["type"][mask_data_xgb_ggf],
        "part":arr_data_xgb_ggf["part"][mask_data_xgb_ggf],
        "lep":arr_data_xgb_ggf["lep"][mask_data_xgb_ggf],
        "mllg_r":arr_data_xgb_ggf["llphoton_refit_m"][mask_data_xgb_ggf],
        "weight":arr_data_xgb_ggf["weight1"][mask_data_xgb_ggf]
        }, index=arr_data_xgb_ggf["index"][mask_data_xgb_ggf])
    
    
    #frames = [df_VBF, df_GGF, df_SMZg, df_DY, df_EWK]
    frames_xgb_ggf = [df_VBF_xgb_ggf, df_GGF_xgb_ggf, df_SMZg_xgb_ggf, df_DY_xgb_ggf, df_EWK_xgb_ggf, df_data_xgb_ggf]
    
    #print(df_VBF)
    #result = pd.concat(frames)
    result_xgb_ggf = pd.concat(frames_xgb_ggf)

#cell ?
#if True:
#    fig = plt.figure(figsize=(12, 10.5))
#    gs = fig.add_gridspec(2, hspace=0, height_ratios=[2,1])
#    (ax0,ax1) = gs.subplots(sharex=True, sharey=False)
#    result = result_xgb_ggf
#    mllg_high = 165
#    mllg_low = 95
#    result_cut_plot = (result.index == 0) & (result.mllg_r <= mllg_high) & (result.mllg_r >= mllg_low)
#    result_cut_plot_data = (result.index == 3) & (result.mllg_r <= mllg_high) & (result.mllg_r >= mllg_low)
#    result_cut_sig = (result.index == 1) & (result.mllg_r <= mllg_high) & (result.mllg_r >= mllg_low)
#    
#    
#    result_cut_plot_ = (result.index == 0) & (result.mllg_r <= 130) & (result.mllg_r >= 120)
#    result_cut_plot_data_ = (result.index == 3) & (result.mllg_r <= 130) & (result.mllg_r >= 120)
#    result_cut_sig_ = (result.index == 1) & (result.mllg_r <= 130) & (result.mllg_r >= 120)
#    
#    plt.rcParams.update({'font.size': 20})
#    val= "_val"
#    nbins = 100
#    range_mask_data = ((arr_data_xgb_ggf["llphoton_refit_m"] > 120) & (arr_data_xgb_ggf["llphoton_refit_m"] < 130)) & mask_data_xgb_ggf 
#    range_mask_GGF = (arr_GGF_xgb_ggf["llphoton_refit_m"] > 120) & (arr_GGF_xgb_ggf["llphoton_refit_m"] < 130) & (arr_GGF_xgb_ggf["weight1"] < 0.6) & mask_GGF_xgb_ggf
#    range_mask_VBF = (arr_VBF_xgb_ggf["llphoton_refit_m"] > 120) & (arr_VBF_xgb_ggf["llphoton_refit_m"] < 130)  & mask_VBF_xgb_ggf
#    range_mask_SMZg = (arr_SMZg_xgb_ggf["llphoton_refit_m"] > 120) & (arr_SMZg_xgb_ggf["llphoton_refit_m"] < 130) & mask_SMZg_xgb_ggf
#    range_mask_DY = (arr_DY_xgb_ggf["llphoton_refit_m"] > 120) & (arr_DY_xgb_ggf["llphoton_refit_m"] < 130) & mask_DY_xgb_ggf
#    range_mask_EWK = (arr_EWK_xgb_ggf["llphoton_refit_m"] > 120) & (arr_EWK_xgb_ggf["llphoton_refit_m"] < 130)  & mask_EWK_xgb_ggf
#    
#    
#    #range_mask_GGF = (arr_GGF_xgb_ggf["weight1"] < 0.6)
#    #range_mask_VBF = (arr_VBF_xgb_ggf["weight1"] < 0.6)
#    #range_mask_SMZg = (arr_SMZg_xgb_ggf["weight1"] < 0.6)
#    #range_mask_DY = (arr_DY_xgb_ggf["weight1"] < 0.6)
#    #range_mask_EWK = (arr_EWK_xgb_ggf["weight1"] < 0.6)
#    
#    
#    ax0.margins(x=0)
#    bkg = sum(result[result_cut_plot]["weight"])
#    data = sum(result[result_cut_plot_data]["weight"])
#    A = data/bkg
#    print(A)
#    
#    S = sum(result[result_cut_plot_data_]["weight"])/sum(result[result_cut_sig_]["weight"])
#    print(S)
#    
#    entries, edges, _ = ax0.hist([arr_SMZg_xgb_ggf["BDT_score"+val][range_mask_SMZg],\
#              arr_DY_xgb_ggf["BDT_score"+val][range_mask_DY], arr_EWK_xgb_ggf["BDT_score"+val][range_mask_EWK]], bins=nbins, \
#             weights = [arr_SMZg_xgb_ggf["weight1"][range_mask_SMZg]*A,\
#                        arr_DY_xgb_ggf["weight1"][range_mask_DY]*A,arr_EWK_xgb_ggf["weight1"][range_mask_EWK]*A],color=["#3f90da","#ffa90e","#92dadd"],\
#             label=['Z+$\\gamma$','Z+FakePhoton stacked','VBSZ+$\\gamma$ stacked'],range=[0,1],alpha=1, stacked=True)
#    
#    lg = Line2D([], [], color='green', lw = 3)    
#    #    h = np.histogram(result[result_cut]['mllg'], bins=nbins, weights =result[result_cut]['weight'])
#    #err_mc = np.sqrt(np.histogram(result[result_cut_plot]['mllg_r'], bins=nbins, weights=(result[result_cut_plot]['weight'])**2)[0])
#    
#    
#    #entries_, edges_, __ = ax0.hist([arr_VBF_xgb_ggf["BDT_score"+val][range_mask_VBF],\
#    #          arr_GGF_xgb_ggf["BDT_score"+val][range_mask_GGF]], bins=nbins,\
#    #         weights = [arr_VBF_xgb_ggf["weight1"][range_mask_VBF]*300,\
#    #                    arr_GGF_xgb_ggf["weight1"][range_mask_GGF]*300], color=['green','red'],label=["VBF","ggF stacked"],histtype=u'step',range=[0,1],alpha=1,linewidth=3, stacked=True)
#    
#    entries_, edges_, __ = ax0.hist(result[result_cut_sig_].ggf_xgb,bins=nbins,weights=result[result_cut_sig_].weight*300,color='red',label="ggF",histtype=u'step',range=[0,1],alpha=1,linewidth=3)
#    if do_vbf or not do_vbf:
#      entries_, edges_, __ = ax0.hist([arr_VBF_xgb_ggf["BDT_score"+val][range_mask_VBF]], bins=nbins,\
#               weights = [arr_VBF_xgb_ggf["weight1"][range_mask_VBF]*300], color=['green'],label=["VBF"],histtype=u'step',range=[0,1],alpha=1,linewidth=3, stacked=True)
#    #entries_, edges_, __ = ax0.hist([arr_GGF_xgb_ggf["BDT_score"+val][range_mask_GGF]], bins=nbins,\
#    #         weights = [arr_GGF_xgb_ggf["weight1"][range_mask_GGF]*300,], color=['red'],label=["ggF"],histtype=u'step',range=[0,1],alpha=1,linewidth=3, stacked=True)
#    bin_edges = np.linspace(0, 1, nbins + 1)
#    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#    bin_widths = bin_edges[1:] - bin_edges[:-1]
#    
#    # --- MC stat uncertainty for stacked background ---
#    smzg_vals = np.asarray(arr_SMZg_xgb_ggf["BDT_score"+val][range_mask_SMZg])
#    dy_vals   = np.asarray(arr_DY_xgb_ggf["BDT_score"+val][range_mask_DY])
#    ewk_vals  = np.asarray(arr_EWK_xgb_ggf["BDT_score"+val][range_mask_EWK])
#    
#    smzg_w = np.asarray(arr_SMZg_xgb_ggf["weight1"][range_mask_SMZg] * A)
#    dy_w   = np.asarray(arr_DY_xgb_ggf["weight1"][range_mask_DY] * A)
#    ewk_w  = np.asarray(arr_EWK_xgb_ggf["weight1"][range_mask_EWK] * A)
#    
#    sumw2_smzg = np.histogram(smzg_vals, bins=edges, weights=smzg_w**2)[0]
#    sumw2_dy   = np.histogram(dy_vals,   bins=edges, weights=dy_w**2)[0]
#    sumw2_ewk  = np.histogram(ewk_vals,  bins=edges, weights=ewk_w**2)[0]
#    
#    err_mc_bkg = np.sqrt(sumw2_smzg + sumw2_dy + sumw2_ewk)
#    
#    
#    vbf_vals = np.asarray(arr_VBF_xgb_ggf["BDT_score"+val][range_mask_VBF])
#    ggf_vals = np.asarray(arr_GGF_xgb_ggf["BDT_score"+val][range_mask_GGF])
#    
#    vbf_w = np.asarray(arr_VBF_xgb_ggf["weight1"][range_mask_VBF])
#    ggf_w = np.asarray(arr_GGF_xgb_ggf["weight1"][range_mask_GGF])
#    
#    sumw2_vbf = np.histogram(vbf_vals, bins=bin_edges, weights=vbf_w**2)[0]
#    sumw2_ggf = np.histogram(ggf_vals, bins=bin_edges, weights=ggf_w**2)[0]
#    err_mc_sig = np.sqrt(sumw2_vbf + sumw2_ggf)
#    #-----------------_#
#    
#    
#    n,bins,patches = ax0.hist(arr_data_xgb_ggf["BDT_score"+val][range_mask_data],\
#                              weights = arr_data_xgb_ggf["weight1"][range_mask_data], color=['green'],histtype=u'step',range=[0,1],alpha=0,linewidth=3, bins=nbins)
#    
#    entries_ = entries_/300
#    err = []
#    for n_ in n:
#        err.append(math.sqrt(n_))
#    err_handle = ax0.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n, yerr=err, fmt="o", color='black', label='Data')
#    
#    
#    mc_total = entries[2] + entries_[1]
#    mc_safe = np.where(mc_total > 0, mc_total, np.nan)
#    err_mc_total = np.sqrt(err_mc_bkg**2 + err_mc_sig**2)
#    
#    
#    ax0.bar(bin_centers,2 * err_mc_bkg,bottom=entries[2] - err_mc_bkg,width=bin_widths,align='center',color='black',alpha=0.12,edgecolor='none',zorder=1)
#    
#    #err = []
#    #for n_ in n:
#    #    err.append(math.sqrt(n_))
#    #ax0.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n, yerr=err, fmt="o", color='black', label='Data')
#    
#    plt.legend()
#    ax0.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True)
#    ax0.tick_params(axis='x', which='minor', bottom=True, labelbottom=True)
#    ax0.tick_params(axis='y', which='minor', left=True, labelleft=True)
#    ax0.xaxis.set_tick_params(labelsize=20)
#    ax0.yaxis.set_tick_params(labelsize=20)
#    #hep.cms.label("Preliminary", loc=1, ax=ax0, data=True, fontsize = 25, rlabel = "138 fb$^{-1}$ (13 TeV) + 62 fb$^{-1}$ (13.6 TeV)")
#    hep.cms.label("", loc=1, ax=ax0, data=True, fontsize = 25, rlabel = "138 fb$^{-1}$ (13 TeV) + 62 fb$^{-1}$ (13.6 TeV)")
#    band_patch=mpatches.Patch(facecolor='black',alpha=0.12,edgecolor='black',label='MC stat. unc.')
#    #legend=ax0.legend(handles=[dot,teal_patch,orange_patch,blue_patch,l,lg],labels = ['Data','VBS Z+$\gamma$','Z+FakePhoton','Z+$\gamma$','Total Signal x 300','VBF x 300', 'Data'], loc="upper right", fontsize=20) #+"+/-"+str(round(A_err*100,2))
#    #legend=ax0.legend(handles=[dot,teal_patch,orange_patch,blue_patch,l],labels = ['Data','VBS Z+$\gamma$','Z+FakePhoton','Z+$\gamma$','Total Signal x 300'], loc="upper right", fontsize=20) #+"+/-"+str(round(A_err*100,2))
#    legend = ax0.legend(handles=[err_handle,teal_patch,orange_patch,blue_patch,l,band_patch],labels = \
#               ['Data','VBS Z+$\\gamma$','Z+Fake photon','Z+$\\gamma$', 'Total signal x 300', 'MC stat. unc.'], loc="upper right", fontsize=20,frameon=False)
#    
#    
#    for text in legend.get_texts():
#        text.set_bbox(dict(facecolor='white', edgecolor='none', pad=5))
#    
#    bdts = [0.91, 0.79, 0.59]
#    ylim = ax0.get_ylim()
#    #bdts = [0.95, 0.86, 0.66] #vbf
#    from matplotlib.pyplot import text
#    #ax0.plot([bdts[0],bdts[0]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
#    #text(bdts[0], 1, "%.2f" % bdts[0], rotation=90, verticalalignment='top')
#    
#    #ax0.text(bdts[0], 1, "%.2f" % bdts[0],transform=ax0.transAxes, fontsize=8, va='top', ha='left')
#    
#    #ax0.plot([bdts[1],bdts[1]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
#    #text(bdts[1], 1, "%.2f" % bdts[1], rotation=90, verticalalignment='top')
#    #ax0.plot([bdts[2],bdts[2]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
#    #text(bdts[2], 1, "%.2f" % bdts[2], rotation=90, verticalalignment='top')
#    ax0.tick_params(axis="x",labelbottom=False)
#    if not do_vbf:
#      ax0.set_ylim(0,2600)
#    else:
#      ax0.set_ylim(0,300)
#    #ax0.set_yscale('log')
#    #ax0.set_ylim(0.001,100000000)
#    
#    fig.add_subplot(ax0)
#    ax1.margins(x=0)
#    ax1.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True)    
#    #ax1.set_ylim(0.5,1.5)
#    
#         #   a1.set_ylim(0.003,2)
#    ax1.axhline(1, linestyle='--',color='black')
#    ax1.axhline(1.5, linestyle='--',color='grey')
#    ax1.axhline(0.5, linestyle='--',color='grey')
#    
#    #ax1.axhline(1, linestyle='--',color='black')
#    #ax1.axhline(1.1, linestyle='--',color='grey')
#    #ax1.axhline(0.9, linestyle='--',color='grey')
#    ax1.set_xlabel("ggF BDT", fontsize = 30)
#    ax1.xaxis.set_tick_params(labelsize=24)
#    ax1.yaxis.set_tick_params(labelsize=24)
#    ax1.tick_params(axis='x', which='minor', bottom=True, labelbottom=True)
#    ax1.tick_params(axis='y', which='minor', left=True, labelleft=True)
#    ax1.yaxis.set_ticks([0.5,1,1.5])
#    ax1.yaxis.set_ticklabels([0.5,1,1.5])
#       # entries_ = [0.0001 if k < 0 else k for k in entries[2]]
#    #yerr_mc = (err_mc*n)/((entries[2])**2)
#    #entries[2][entries[2]<0]=0.00001
#    #a1.
#    #a1.errorbar(edges[:-1]+ 0.5*(edges[1:] - edges[:-1]), entries[0]/entries_[0], yerr=1, fmt=".", color='darkred')
#    #entries[2][entries[2]<0]=0.00001
#    #ax1.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n/ (entries[2]+entries_[1]), yerr=err/(entries[2]+entries_[1]), fmt="o", color='black')
#    err = np.asarray(err)
#    ratio = n / mc_safe
#    ratio_err = err / mc_safe
#    mc_band = err_mc_total / mc_safe
#    
#    #ax1.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n/mc_total, yerr=err/mc_total, fmt="o", color='black')
#    ax1.errorbar(bin_centers,ratio,yerr=ratio_err,fmt="o",color="black",zorder=2)
#    ax1.bar(bin_centers,2 * mc_band,bottom=1 - mc_band,width=bin_widths,align='center',color='black',alpha=0.12,edgecolor='none',zorder=0)
#     
#    #ax1.errorbar(edges[:-1]+ 0.5*(edges[1:] - edges[:-1]), entries[1]/entries_[2], yerr=0, fmt="-", color='red')
#    #ax1.fill_between(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), 1-yerr_mc, 1+yerr_mc, alpha=0.1, edgecolor='none', facecolor='black',
#    #linewidth=4,step='pre')
#    #    ax.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), np.ones(len(n)), yerr=(err_mc*n)/((entries[2])**2), color='grey',fillstye='full')
#    ax1.set_ylabel('Data/MC', fontsize = 20)    
#    ax1.set_ylim(0.02,2.05)
#    #box = a1.get_position()
#    #a1.set_position([box.x0, box.y0, box.width, box.height])
#    
#    #plt.ylim([0,400])
#    #plt.ylim([0.01,1000])
#    ylim = plt.ylim()
#    
#    #from matplotlib.pyplot import text
#    #ax.set_yscale('log')
#    #bdts = [0.85,0.64, 0.42] #ggf
#    #bdts = [0.91, 0.82, 0.61]
#    #bdts = [0.91, 0.79, 0.59]
#    #bdts = [0.95, 0.86, 0.66] #vbf
#    bdts = [0.94,0.83,0.57]
#    plt.plot([bdts[0],bdts[0]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
#    #text(bdts[0], 1, "%.2f" % bdts[0], rotation=90, verticalalignment='center')
#    plt.plot([bdts[1],bdts[1]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
#    #text(bdts[1], 1, "%.2f" % bdts[1], rotation=90, verticalalignment='center')
#    plt.plot([bdts[2],bdts[2]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
#    #text(bdts[2], 1, "%.2f" % bdts[2], rotation=90, verticalalignment='center')
#    bdt_edges = [0.0, 0.57, 0.83, 0.94, 1.0]
#    text_yoffsets = [1750, 1240, 1000, 1000]
#    for ilabel in range(4):
#      x_pos = (bdt_edges[ilabel]+bdt_edges[ilabel+1])/2.0
#      plt.text(x_pos, text_yoffsets[ilabel], f'ggF{4-ilabel}', fontsize=18.0, 
#               ha='center')
#    
#    plt.ylabel("Events / 0.01",fontsize=20)
#    
#    #plt.savefig("BDT_ggF_unblind.png")
#    plt.savefig("BDT_score_ggF_unblind.pdf")

#VBF
if True:
    fig = plt.figure(figsize=(12, 10.5))
    gs = fig.add_gridspec(2, hspace=0, height_ratios=[2,1])
    (ax0,ax1) = gs.subplots(sharex=True, sharey=False)
    result = result_xgb_ggf
    mllg_high = 165
    mllg_low = 95
    result_cut_plot = (result.index == 0) & (result.mllg_r <= mllg_high) & (result.mllg_r >= mllg_low)
    result_cut_plot_data = (result.index == 3) & (result.mllg_r <= mllg_high) & (result.mllg_r >= mllg_low)
    result_cut_sig = (result.index == 1) & (result.mllg_r <= mllg_high) & (result.mllg_r >= mllg_low)
    
    
    result_cut_plot_ = (result.index == 0) & (result.mllg_r <= 130) & (result.mllg_r >= 120)
    result_cut_plot_data_ = (result.index == 3) & (result.mllg_r <= 130) & (result.mllg_r >= 120)
    result_cut_sig_ = (result.index == 1) & (result.mllg_r <= 130) & (result.mllg_r >= 120)
    
    plt.rcParams.update({'font.size': 20})
    val= "_val"
    nbins = 100
    range_mask_data = ((arr_data_xgb_ggf["llphoton_refit_m"] > 120) & (arr_data_xgb_ggf["llphoton_refit_m"] < 130)) & mask_data_xgb_ggf 
    range_mask_GGF = (arr_GGF_xgb_ggf["llphoton_refit_m"] > 120) & (arr_GGF_xgb_ggf["llphoton_refit_m"] < 130) & (arr_GGF_xgb_ggf["weight1"] < 0.6) & mask_GGF_xgb_ggf
    range_mask_VBF = (arr_VBF_xgb_ggf["llphoton_refit_m"] > 120) & (arr_VBF_xgb_ggf["llphoton_refit_m"] < 130)  & mask_VBF_xgb_ggf
    range_mask_SMZg = (arr_SMZg_xgb_ggf["llphoton_refit_m"] > 120) & (arr_SMZg_xgb_ggf["llphoton_refit_m"] < 130) & mask_SMZg_xgb_ggf
    range_mask_DY = (arr_DY_xgb_ggf["llphoton_refit_m"] > 120) & (arr_DY_xgb_ggf["llphoton_refit_m"] < 130) & mask_DY_xgb_ggf
    range_mask_EWK = (arr_EWK_xgb_ggf["llphoton_refit_m"] > 120) & (arr_EWK_xgb_ggf["llphoton_refit_m"] < 130)  & mask_EWK_xgb_ggf
    
    
    #range_mask_GGF = (arr_GGF_xgb_ggf["weight1"] < 0.6)
    #range_mask_VBF = (arr_VBF_xgb_ggf["weight1"] < 0.6)
    #range_mask_SMZg = (arr_SMZg_xgb_ggf["weight1"] < 0.6)
    #range_mask_DY = (arr_DY_xgb_ggf["weight1"] < 0.6)
    #range_mask_EWK = (arr_EWK_xgb_ggf["weight1"] < 0.6)
    
    
    ax0.margins(x=0)
    bkg = sum(result[result_cut_plot]["weight"])
    data = sum(result[result_cut_plot_data]["weight"])
    A = data/bkg
    print(A)
    
    S = sum(result[result_cut_plot_data_]["weight"])/sum(result[result_cut_sig_]["weight"])
    print(S)
    
    #plot background stack
    entries, edges, _ = ax0.hist([arr_SMZg_xgb_ggf["BDT_score"+val][range_mask_SMZg],\
              arr_DY_xgb_ggf["BDT_score"+val][range_mask_DY], arr_EWK_xgb_ggf["BDT_score"+val][range_mask_EWK]], bins=nbins, \
             weights = [arr_SMZg_xgb_ggf["weight1"][range_mask_SMZg]*A,\
                        arr_DY_xgb_ggf["weight1"][range_mask_DY]*A,arr_EWK_xgb_ggf["weight1"][range_mask_EWK]*A],color=["#3f90da","#ffa90e","#92dadd"],\
             label=['Z+$\\gamma$','Z+FakePhoton stacked','VBSZ+$\\gamma$ stacked'],range=[0,1],alpha=1, stacked=True)
    
    lg = Line2D([], [], color='green', lw = 3)    
    #    h = np.histogram(result[result_cut]['mllg'], bins=nbins, weights =result[result_cut]['weight'])
    #err_mc = np.sqrt(np.histogram(result[result_cut_plot]['mllg_r'], bins=nbins, weights=(result[result_cut_plot]['weight'])**2)[0])
    
    #plot signal
    signal_sf = 250
    if do_vbf:
        signal_sf = 130
    
    #entries_, edges_, __ = ax0.hist([arr_VBF_xgb_ggf["BDT_score"+val][range_mask_VBF],\
    #          arr_GGF_xgb_ggf["BDT_score"+val][range_mask_GGF]], bins=nbins,\
    #         weights = [arr_VBF_xgb_ggf["weight1"][range_mask_VBF]*signal_sf,\
    #                    arr_GGF_xgb_ggf["weight1"][range_mask_GGF]*signal_sf], color=['green','red'],label=["VBF","ggF stacked"],histtype=u'step',range=[0,1],alpha=1,linewidth=3, stacked=True)
    
    entries_, edges_, __ = ax0.hist(result[result_cut_sig_].ggf_xgb,bins=nbins,weights=result[result_cut_sig_].weight*signal_sf,color='red',label="ggF",histtype=u'step',range=[0,1],alpha=1,linewidth=3)
    if do_vbf or not do_vbf:
        entries_, edges_, __ = ax0.hist([arr_VBF_xgb_ggf["BDT_score"+val][range_mask_VBF]], bins=nbins,\
                 weights = [arr_VBF_xgb_ggf["weight1"][range_mask_VBF]*signal_sf], color=['green'],label=["VBF"],histtype=u'step',range=[0,1],alpha=1,linewidth=3, stacked=True)
    #entries_, edges_, __ = ax0.hist([arr_GGF_xgb_ggf["BDT_score"+val][range_mask_GGF]], bins=nbins,\
    #         weights = [arr_GGF_xgb_ggf["weight1"][range_mask_GGF]*300,], color=['red'],label=["ggF"],histtype=u'step',range=[0,1],alpha=1,linewidth=3, stacked=True)
    bin_edges = np.linspace(0, 1, nbins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    
    # --- MC stat uncertainty for stacked background ---
    smzg_vals = np.asarray(arr_SMZg_xgb_ggf["BDT_score"+val][range_mask_SMZg])
    dy_vals   = np.asarray(arr_DY_xgb_ggf["BDT_score"+val][range_mask_DY])
    ewk_vals  = np.asarray(arr_EWK_xgb_ggf["BDT_score"+val][range_mask_EWK])
    
    smzg_w = np.asarray(arr_SMZg_xgb_ggf["weight1"][range_mask_SMZg] * A)
    dy_w   = np.asarray(arr_DY_xgb_ggf["weight1"][range_mask_DY] * A)
    ewk_w  = np.asarray(arr_EWK_xgb_ggf["weight1"][range_mask_EWK] * A)
    
    sumw2_smzg = np.histogram(smzg_vals, bins=edges, weights=smzg_w**2)[0]
    sumw2_dy   = np.histogram(dy_vals,   bins=edges, weights=dy_w**2)[0]
    sumw2_ewk  = np.histogram(ewk_vals,  bins=edges, weights=ewk_w**2)[0]
    
    err_mc_bkg = np.sqrt(sumw2_smzg + sumw2_dy + sumw2_ewk)
    
    
    vbf_vals = np.asarray(arr_VBF_xgb_ggf["BDT_score"+val][range_mask_VBF])
    ggf_vals = np.asarray(arr_GGF_xgb_ggf["BDT_score"+val][range_mask_GGF])
    
    vbf_w = np.asarray(arr_VBF_xgb_ggf["weight1"][range_mask_VBF])
    ggf_w = np.asarray(arr_GGF_xgb_ggf["weight1"][range_mask_GGF])
    
    sumw2_vbf = np.histogram(vbf_vals, bins=bin_edges, weights=vbf_w**2)[0]
    sumw2_ggf = np.histogram(ggf_vals, bins=bin_edges, weights=ggf_w**2)[0]
    err_mc_sig = np.sqrt(sumw2_vbf + sumw2_ggf)
    #-----------------_#
    
    
    n,bins,patches = ax0.hist(arr_data_xgb_ggf["BDT_score"+val][range_mask_data],\
                              weights = arr_data_xgb_ggf["weight1"][range_mask_data], color=['green'],histtype=u'step',range=[0,1],alpha=0,linewidth=3, bins=nbins)
    
    entries_ = entries_/300
    err = []
    for n_ in n:
        err.append(math.sqrt(n_))
    err_handle = ax0.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n, yerr=err, fmt="o", color='black', label='Data')
    
    
    mc_total = entries[2] + entries_[1]
    mc_safe = np.where(mc_total > 0, mc_total, np.nan)
    err_mc_total = np.sqrt(err_mc_bkg**2 + err_mc_sig**2)
    
    
    ax0.bar(bin_centers,2 * err_mc_bkg,bottom=entries[2] - err_mc_bkg,width=bin_widths,align='center',color='black',alpha=0.12,edgecolor='none',zorder=1)
    
    #err = []
    #for n_ in n:
    #    err.append(math.sqrt(n_))
    #ax0.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n, yerr=err, fmt="o", color='black', label='Data')
    
    plt.legend()
    ax0.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True)
    ax0.tick_params(axis='x', which='minor', bottom=True, labelbottom=True)
    ax0.tick_params(axis='y', which='minor', left=True, labelleft=True)
    ax0.xaxis.set_tick_params(labelsize=20)
    ax0.yaxis.set_tick_params(labelsize=20)
    #hep.cms.label("Preliminary", loc=1, ax=ax0, data=True, fontsize = 25, rlabel = "138 fb$^{-1}$ (13 TeV) + 62 fb$^{-1}$ (13.6 TeV)")
    hep.cms.label("", loc=1, ax=ax0, data=True, fontsize = 25, rlabel = "138 fb$^{-1}$ (13 TeV) + 62 fb$^{-1}$ (13.6 TeV)")
    band_patch=mpatches.Patch(facecolor='black',alpha=0.12,edgecolor='black',label='MC stat. unc.')
    #legend=ax0.legend(handles=[dot,teal_patch,orange_patch,blue_patch,l,lg],labels = ['Data','VBS Z+$\gamma$','Z+FakePhoton','Z+$\gamma$','Total Signal x 300','VBF x 300', 'Data'], loc="upper right", fontsize=20) #+"+/-"+str(round(A_err*100,2))
    #legend=ax0.legend(handles=[dot,teal_patch,orange_patch,blue_patch,l],labels = ['Data','VBS Z+$\gamma$','Z+FakePhoton','Z+$\gamma$','Total Signal x 300'], loc="upper right", fontsize=20) #+"+/-"+str(round(A_err*100,2))

    #if not do_vbf:
    #    legend = ax0.legend(handles=[err_handle,teal_patch,orange_patch,blue_patch,l,band_patch],labels = \
    #               ['Data','VBS Z+$\\gamma$','Z+Fake photon','Z+$\\gamma$', 'Total signal x 300', 'MC stat. unc.'], loc="upper right", fontsize=20,frameon=False)
    if do_vbf or not do_vbf:
        legend = ax0.legend(handles=[err_handle,teal_patch,orange_patch,blue_patch,l,lg,band_patch],labels = \
                   ['Data','VBS Z+$\\gamma$','Z+Fake photon','Z+$\\gamma$', f'Total signal x {signal_sf}', f'VBF x {signal_sf}', 'MC stat. unc.'], loc="upper right", fontsize=20,frameon=False)
    
    
    for text in legend.get_texts():
        text.set_bbox(dict(facecolor='white', edgecolor='none', pad=5))
    
    bdts = [0.91, 0.79, 0.59]
    ylim = ax0.get_ylim()
    #bdts = [0.95, 0.86, 0.66] #vbf
    from matplotlib.pyplot import text
    #ax0.plot([bdts[0],bdts[0]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
    #text(bdts[0], 1, "%.2f" % bdts[0], rotation=90, verticalalignment='top')
    
    #ax0.text(bdts[0], 1, "%.2f" % bdts[0],transform=ax0.transAxes, fontsize=8, va='top', ha='left')
    
    #ax0.plot([bdts[1],bdts[1]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
    #text(bdts[1], 1, "%.2f" % bdts[1], rotation=90, verticalalignment='top')
    #ax0.plot([bdts[2],bdts[2]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
    #text(bdts[2], 1, "%.2f" % bdts[2], rotation=90, verticalalignment='top')
    ax0.tick_params(axis="x",labelbottom=False)
    if not do_vbf:
      ax0.set_ylim(0,2600)
    else:
      ax0.set_ylim(0,300)
    #ax0.set_yscale('log')
    #ax0.set_ylim(0.001,100000000)
    
    fig.add_subplot(ax0)
    ax1.margins(x=0)
    ax1.tick_params(bottom=True, labelbottom=True, left=True, labelleft=True)    
    #ax1.set_ylim(0.5,1.5)
    
         #   a1.set_ylim(0.003,2)
    ax1.axhline(1, linestyle='--',color='black')
    ax1.axhline(1.5, linestyle='--',color='grey')
    ax1.axhline(0.5, linestyle='--',color='grey')
    
    #ax1.axhline(1, linestyle='--',color='black')
    #ax1.axhline(1.1, linestyle='--',color='grey')
    #ax1.axhline(0.9, linestyle='--',color='grey')
    if not do_vbf:
        ax1.set_xlabel("ggF BDT score", fontsize = 30)
    else:
        ax1.set_xlabel("VBF BDT score", fontsize = 30)
    ax1.xaxis.set_tick_params(labelsize=24)
    ax1.yaxis.set_tick_params(labelsize=24)
    ax1.tick_params(axis='x', which='minor', bottom=True, labelbottom=True)
    ax1.tick_params(axis='y', which='minor', left=True, labelleft=True)
    ax1.yaxis.set_ticks([0.5,1,1.5])
    ax1.yaxis.set_ticklabels([0.5,1,1.5])
       # entries_ = [0.0001 if k < 0 else k for k in entries[2]]
    #yerr_mc = (err_mc*n)/((entries[2])**2)
    #entries[2][entries[2]<0]=0.00001
    #a1.
    #a1.errorbar(edges[:-1]+ 0.5*(edges[1:] - edges[:-1]), entries[0]/entries_[0], yerr=1, fmt=".", color='darkred')
    #entries[2][entries[2]<0]=0.00001
    #ax1.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n/ (entries[2]+entries_[1]), yerr=err/(entries[2]+entries_[1]), fmt="o", color='black')
    err = np.asarray(err)
    ratio = n / mc_safe
    ratio_err = err / mc_safe
    mc_band = err_mc_total / mc_safe
    
    #ax1.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), n/mc_total, yerr=err/mc_total, fmt="o", color='black')
    ax1.errorbar(bin_centers,ratio,yerr=ratio_err,fmt="o",color="black",zorder=2)
    ax1.bar(bin_centers,2 * mc_band,bottom=1 - mc_band,width=bin_widths,align='center',color='black',alpha=0.12,edgecolor='none',zorder=0)
     
    #ax1.errorbar(edges[:-1]+ 0.5*(edges[1:] - edges[:-1]), entries[1]/entries_[2], yerr=0, fmt="-", color='red')
    #ax1.fill_between(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), 1-yerr_mc, 1+yerr_mc, alpha=0.1, edgecolor='none', facecolor='black',
    #linewidth=4,step='pre')
    #    ax.errorbar(bins[:-1]+ 0.5*(bins[1:] - bins[:-1]), np.ones(len(n)), yerr=(err_mc*n)/((entries[2])**2), color='grey',fillstye='full')
    ax1.set_ylabel('Data/MC', fontsize = 20)    
    ax1.set_ylim(0.02,2.05)
    #box = a1.get_position()
    #a1.set_position([box.x0, box.y0, box.width, box.height])
    
    #plt.ylim([0,400])
    #plt.ylim([0.01,1000])
    ylim = plt.ylim()
    
    #from matplotlib.pyplot import text
    #ax.set_yscale('log')
    #bdts = [0.85,0.64, 0.42] #ggf
    #bdts = [0.91, 0.82, 0.61]
    #bdts = [0.91, 0.79, 0.59]
    #bdts = [0.95, 0.86, 0.66] #vbf
    bdts = [0.94,0.83,0.57]
    if do_vbf:
        bdts = [0.91,0.81,0.48]
    plt.plot([bdts[0],bdts[0]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
    #text(bdts[0], 1, "%.2f" % bdts[0], rotation=90, verticalalignment='center')
    plt.plot([bdts[1],bdts[1]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
    #text(bdts[1], 1, "%.2f" % bdts[1], rotation=90, verticalalignment='center')
    plt.plot([bdts[2],bdts[2]],[ylim[0],ylim[1]], c='black',ls='--',alpha=0.5)
    #text(bdts[2], 1, "%.2f" % bdts[2], rotation=90, verticalalignment='center')
    if not do_vbf:
        text_xoffsets = [0.285, 0.74, 0.885, 0.97]
        text_yoffsets = [1750, 900, 900, 900]
        for ilabel in range(4):
            plt.text(text_xoffsets[ilabel], text_yoffsets[ilabel], f'ggF{4-ilabel}', 
                     fontsize=18.0, ha='center')
    else:
        text_xoffset = [0.24, 0.645, 0.86, 0.955]
        text_yoffsets = [220, 105, 105, 105]
        for ilabel in range(4):
            plt.text(text_xoffset[ilabel], text_yoffsets[ilabel], 
                     f'VBF{4-ilabel}', 
                     fontsize=18.0, ha='center')
    
    plt.ylabel("Events / 0.01",fontsize=20)
    
    #plt.savefig("BDT_ggF_unblind.png")
    if not do_vbf:
        plt.savefig("BDT_score_ggF_unblind.pdf")
    else:
        plt.savefig("BDT_score_VBF_unblind.pdf")

    #generate text file for hepdata
    sumw_smzg = np.histogram(smzg_vals, bins=edges, weights=smzg_w)[0]
    sumw_dy   = np.histogram(dy_vals,   bins=edges, weights=dy_w)[0]
    sumw_ewk  = np.histogram(ewk_vals,  bins=edges, weights=ewk_w)[0]
    sumw_vbf = np.histogram(vbf_vals, bins=bin_edges, weights=vbf_w)[0]
    sumw_ggf = np.histogram(ggf_vals, bins=bin_edges, weights=ggf_w)[0]
    hepdata_data = concatenate_flat([
        bin_centers, sumw_smzg, sumw2_smzg, sumw_dy, sumw2_dy, 
        sumw_ewk, sumw2_ewk, sumw_vbf, sumw2_vbf, sumw_ggf, sumw2_ggf, 
        sumw_smzg+sumw_dy+sumw_ewk, err_mc_bkg**2, sumw_ggf+sumw_vbf, 
        err_mc_sig**2, n, err, ratio, ratio_err, mc_band])
    if not do_vbf:
      np.savetxt("BDT_score_ggF_unblind.txt", hepdata_data)
    else:
      np.savetxt("BDT_score_VBF_unblind.txt", hepdata_data)

