# Get dictionaries
#dx extract_dataset project-GvFxJ08J95KXx97XFz8g2X2g:record-Gzf082QJB0k8bF6ByK1xzfb1 -ddd --delimiter ","

dict="ref/Cohort_Browser.data_dictionary.csv"
export PATH="$HOME/.local/bin:$PATH"

echo "eid" > field_name_sajanthDYT_PD.txt

#Added death registry (40001/2), PRS (26260/1), Date & Source of PD report (42032/3), Assessment Centre Visit Date (53) to standard set
fields=("20002" "41270" "20111" "20110" "20107" "22001" "22006" "22009" "22000" "34" "21022" "22021" "22019" "22027" "22189" "40001" "40002" "26260" "26261" "42032" "42033" "53") 

for id in "${fields[@]}"; do
    awk -F',' '{print $2}' $dict | grep -E "^p$id([^0-9]|$)"  >> field_name_sajanthDYT_PD.txt;
done

dx upload field_name_sajanthDYT_PD.txt --destination tabular_data/

# Then run Table Exporter 