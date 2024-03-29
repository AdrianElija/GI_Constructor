#Capstone Neighborhood Change Index Construction Script Tool
### imports
import numpy as np
import pandas as pd
import arcpy
import sys
import os
# Get input parameters
Parcels_folder = arcpy.GetParameterAsText(0)
CENACS_folder = arcpy.GetParameterAsText(1)
Geography_Baseline = arcpy.GetParameterAsText(2)
Study_subject = arcpy.GetParameterAsText(3)
Output_File_Location = arcpy.GetParameterAsText(4)
# Get list of feature classes in the input Geodatabase from Parcels_folder
arcpy.env.workspace = Parcels_folder
Parcels_list = arcpy.ListFeatureClasses()
arcpy.AddMessage(f"Found {len(Parcels_list)} Parcels files:")
arcpy.AddMessage(Parcels_list)
if len(Parcels_list) == 0:
    arcpy.AddMessage("No feature classes found in workspace")
else:
    arcpy.AddMessage("Found feature classes:")
    for fc in Parcels_list:
        arcpy.AddMessage(fc)

# Find all Polygon feature classes in all File Geodatabases in CENACS_folder
arcpy.env.workspace = CENACS_folder
CENACS_list = arcpy.ListFeatureClasses()

# Print the list of files for debugging purposes
arcpy.AddMessage(f"Found {len(CENACS_list)} CENACS files:")
arcpy.AddMessage(CENACS_list)

if len(CENACS_list) == 0:
    arcpy.AddMessage("No feature classes found in workspace")
else:
    arcpy.AddMessage("Found feature classes:")
    for fc in CENACS_list:
        arcpy.AddMessage(fc)
arcpy.AddMessage(CENACS_list)
arcpy.AddMessage(Parcels_list)

def createanalysisgeos(CENACS_list, Output_File_Location,Study_subject):
    studybuffer = os.path.join(Output_File_Location, "\study_buffer")
    buffer_distance_or_field = "0.5 Miles"
    sideType = "FULL"
    endType = "ROUND"
    #Create 0.5 mile buffer
    arcpy.Buffer_analysis(Study_subject, studybuffer, buffer_distance_or_field, sideType, endType)
    for cenacs in CENACS_list:
        # Select by location to select block groups within geography
        arcpy.SelectLayerByLocation_management(cenacs, "INTERSECT", Geography_Baseline)
        # Save features selected by Geography_Baseline to output location
        out_feature_class1 = os.path.join(Output_File_Location, f"\{cenacs}_baseline")
        arcpy.CopyFeatures_management(cenacs, out_feature_class1)
        # Select by Location: intersect 0.5 mile buffer:
        arcpy.SelectLayerByLocation_management(cenacs, "INTERSECT", studybuffer)
        # Save features selected by 0.5 mile buffer to output location
        studycenacs = os.path.join(Output_File_Location, f"\{cenacs}_study_buffer")
        arcpy.CopyFeatures_management(cenacs, studycenacs)

#Prep the CENACS

def modcenacs(CENACS_list, Output_File_Location):
    arcpy.AddMessage(arcpy.env.workspace)
    # Define field names and data types
    field1 = "PCT_MIN"
    field2 = "PCT_ED_COLLEGE"
    field_type = "DOUBLE"
    for cenacs in CENACS_list:
        # Add fields
        arcpy.AddField_management(cenacs, field1, field_type)
        arcpy.AddField_management(cenacs, field2, field_type)
        # Calculate fields
        expression1 = "((!TOTALPOP!-!WHITE!)/!TOTALPOP!)*100.0"
        arcpy.CalculateField_management(cenacs, field1, expression1, "PYTHON_9.3")
        expression2 = "!ED_COLLEGE!/!TOTALPOP!"
        arcpy.CalculateField_management(cenacs, field2, expression2, "PYTHON_9.3")
        # Save modified dataset to output location
        arcpy.env.overwriteOutput = True

def modparcels(Parcels_list, Output_File_Location, CENACS_list):
    arcpy.AddMessage(arcpy.env.workspace)
    # Define where clauses to select residential and nonresidential parcels
    where_clauses = ["LU_RES = 'YES'", "NOT LU_RES = 'YES'"]
    out_tables = []
    for parcels in Parcels_list:
        arcpy.AddMessage(parcels)
        arcpy.SelectLayerByLocation_management(parcels, 'HAVE_THEIR_CENTER_IN', Geography_Baseline)
        # Loop over the where clauses to select residential and nonresidential parcels
        summary_tables = []
        for j, where_clause in enumerate(where_clauses):
            arcpy.SelectLayerByAttribute_management(parcels, "NEW_SELECTION", where_clause)
            # Save selected parcels to output location
            out_feature_class = os.path.join(Output_File_Location, f"\\{os.path.basename(parcels)}_{'residential' if j == 0 else 'nonresidential'}out_test")
            arcpy.CopyFeatures_management(parcels, out_feature_class)
            # Summarize parcels attributes within census block groups
            for CENACS in CENACS_list:
                in_polygons = CENACS
                in_sum_features = out_feature_class  # Use the selected parcel layer
                out_feature_class2 = os.path.join(Output_File_Location, f"{os.path.basename(parcels)}_{'residential' if j == 0 else 'nonresidential'}_summary")
                out_feature_class2 = arcpy.ValidateTableName(out_feature_class2, Output_File_Location)
                out_feature_class2 = arcpy.CreateUniqueName(out_feature_class2, Output_File_Location)
                out_feature_class2 = arcpy.CreateFeatureclass_management(Output_File_Location,
                                                                         os.path.basename(out_feature_class2),
                                                                         "POLYGON", spatial_reference=CENACS,
                                                                         has_m="DISABLED", has_z="DISABLED")
                sum_fields = ["JV"]
                arcpy.SummarizeWithin_analysis(in_polygons, in_sum_features, out_feature_class2, "KEEP_ALL", f"{sum_fields[0]} SUM;{sum_fields[0]} MEAN")
                # Append the output summary table to the list of summary tables for the current year
                summary_tables.append(out_feature_class2)
                # Append the summary tables for the current year to the list of output tables
            out_tables.append(summary_tables)
            # Print a message for each processed feature class
            arcpy.AddMessage(f"Processed {os.path.basename(parcels)})")
            # Flatten the list of output tables
        out_tables = [table for year in out_tables for table in year]
        # Return the list of output tables
        return out_tables


def datacalcandpctchg(CENACS_list, Output_File_Location):
    output_tables = []
    for CENACS in CENACS_list:
        ##Summary Statistics
        out_table = os.path.join(Output_File_Location, f"{os.path.basename(CENACS)}indicators_summary")
        # On Modded CENACS that now has the summarized JV values for each from parcels
        arcpy.Statistics_analysis(CENACS, out_table, [["MEDHHINC", "MEAN"], ["Res_JV", "MEAN"], ["NonRes_JV", "MEAN"], ["PCT_MIN", "MEAN"], ["PCT_ED_COLLEGE", "MEAN"]])

        ## Calc percent change
        #Over 10 years (CENACS: 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
        df = pd.read_table(out_table, delimiter=",")

        # Calculate the percent change for each indicator
        pct_chg_cols = []
        for year in range(11, 21):
            col_name = "JV{}".format(year)
            pct_chg_col_name = "pct_chg_{}".format(year)
            df[pct_chg_col_name] = (df[col_name] - df["JV10"]) / df["JV10"] * 100
            pct_chg_cols.append(pct_chg_col_name)

        # Summarize the percent change for each indicator
        summary_cols = []
        for col_name in pct_chg_cols:
            summary_col_name = "sum_{}".format(col_name)
            df[summary_col_name] = df[col_name].sum()
            summary_cols.append(summary_col_name)

        # Append the summary table to the output_tables list
        case_field = "ID"
        summary_table = df[["OBJECTID", case_field] + summary_cols]
        output_tables.append(summary_table)

    # Concatenate all the output tables into a single table
    Total_output_table = pd.concat(output_tables)
    df = Total_output_table
    return df

## Index Construction
def indexconstructor(indicators_table, Output_File_Location):
    def indexconstructor(indicators_table, Output_File_Location):
        # Define the function to calculate the weighted score
        def calculate_weighted_score(row):
            pct_min_weight = 0.4
            pct_edcollege_weight = 0.2
            hhmedinc_weight = 0.2
            res_jv_weight = 0.1
            nonres_jv_weight = 0.1
            gi_score = row['GI_Score']
            pct_min = row['Percent Minority Population']
            pct_edcollege = row['Percent Population with Bachelor\'s Degree or Higher']
            hhmedinc = row['Average Median Income']
            res_jv = row['res_JV']
            nonres_jv = row['nonres_JV']
            weighted_score = (gi_score * (
                        pct_min_weight * pct_min + pct_edcollege_weight * pct_edcollege + hhmedinc_weight * hhmedinc + res_jv_weight * res_jv + nonres_jv_weight * nonres_jv)) / 100
            return weighted_score

        # Add field to the indicators_table for the weighted score
        indicators_table['GI_Weighted_Score'] = indicators_table.apply(calculate_weighted_score, axis=1)

        # Normalize the data using min-max normalization
        indicators = ['Percent Minority Population', 'Average Median Income', 'Average Property Values',
                      'Educational Attainment']
        for indicator in indicators:
            indicators_table[indicator] = (indicators_table[indicator] - indicators_table[indicator].min()) / (
                        indicators_table[indicator].max() - indicators_table[indicator].min())

        # Save the normalized table to a CSV file
        normalized_table = os.path.join(Output_File_Location, 'normalized_table.csv')
        indicators_table.to_csv(normalized_table, index=False)


######## Calculate the weighted value ###########

######## Sum total weighted scores ############

######## Main #########
# Create the file geodatabase at the output location
#arcpy.CreateFileGDB_management(Output_File_Location, "output.gdb")

modcenacs(CENACS_list, Output_File_Location)
createanalysisgeos(CENACS_list, Output_File_Location, Study_subject)
# Run the modparcels function to summarize parcel data by census block group
parcel_summaries = modparcels(Parcels_list, Output_File_Location, CENACS_list)

# Run the datacalcandpctchg function to calculate data changes and percentages
datacalcandpctchg(parcel_summaries)
indicators_table = datacalcandpctchg(Output_File_Location, CENACS_list)
indexconstructor(indicators_table, Output_File_Location)
