#This script is to make the control region plots for all categories. 
#The settings will make all of the plots, not just sideband
#To change this, change the setting "all" in the line with the executable file
#IT CAN BE ANYTHING BUT BLANK

#Categories are 0 = ttH leptonic, 1 = ttH hadronic, 2 = ZH MET, 3 = WH/ZH 3l, 4 = VBF, 5 = ggF
categories=(0 1 2 3 4 5)

#Making sure the correct environment is loaded
source set_env.sh

#Loops through each category and runs the control region plots consecutively NOT CONCURRENTLY
#This is done to prevent a lot of competition over disk read on cms11
for cat in ${categories[@]}
do
  nohup run/zgamma/plot_an_controlregions.exe all $cat > output/output_TP_$cat.txt &
done


