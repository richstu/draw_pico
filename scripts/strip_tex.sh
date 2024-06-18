#category=("ttH_lep" "ttH_had" "ZH_MET" "WH_3l" "VBF" "ggF")
#flavor=("ee" "mumu" "ll")
#"ll")

export txt_out_file=$1
echo $txt_out_file

txt_out_file=${txt_out_file%.tex}.txt
echo $txt_out_file

#mkdir ./txt_files
mkdir output
cp tables/$1 ./output/$txt_out_file 

vi -c ":\$-0,\$d" -c ":1,12d" -c ":wq" ./output/$txt_out_file
vi -c ":%s/preview/table/g" -c ":wq" ./output/$txt_out_file
