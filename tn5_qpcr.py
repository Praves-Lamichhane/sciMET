import math
import numpy as np
import re
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy.stats
from sklearn.linear_model import LinearRegression


########################################
## READ THE DATA FROM THE CSV file
#######################################

df = pd.read_csv("2020-03-03_173902.csv", sep=",") ## can add "index_col="Name"" to set default index just like "set_index="""
df = df.rename(columns={"Sample Name":"Samples"})
df_qubit = pd.read_csv("transposome_production_n48.csv", sep=",")
df_qubit = df_qubit[["Samples", "dsDNA (ng/uL)"]]
# pd.set_option("display.max_rows", 300)
pd.set_option("display.max_rows", 400) # this sets up the number of rows displayed
pd.set_option("display.float_format", lambda x: '%.3f' % x)

#########################################
## CLEAN THE DATA--REMOVE THE COLUMNS NOT NEEDED
###########################################

df = df.drop(columns=['Omit', 'Well Position', 'Target Name', 'Ct Mean', 'Task',
       'Reporter', 'Quencher', 'Quantity', 'Ct SD', 'dsDNA (ng/uL)', 'Well',
       'Quantity Mean', 'Quantity SD', 'Y-Intercept', 'R(superscript 2)',
       'Slope', 'Efficiency', 'Automatic Ct Threshold', 'Ct Threshold',
       'Automatic Baseline', 'Baseline Start', 'Baseline End', 'Comments',
       'Tm1', 'Tm2', 'Tm3']) ## clean the columns

#################################
## Calculate the mean and standard deviation
##################################
df['CT'] = pd.to_numeric(df['CT'], errors='coerce') ## convert object values to float

mean_CT = df["CT"].mean()
std_CT = df["CT"].std()

#################################
## Filter the rows by Dixon's Q test (not implemented)
##################################

# df = df.sort_values(by="CT")

filt = df["CT"] < (mean_CT + 3*std_CT)
df = df[filt]
# filt1 = df["CT"] < 10
# df = df[filt1]

##################################
## Calculating the mean of each Sample
##################################

new_df = df.groupby("Samples").mean() ## this indexes the Sample Name column
# print(new_df.index) ## checking if the triplicates had resolved
## notice by this time the negative control got filtered out
new_df = new_df.reset_index()
# new_df = new_df.set_index([pd.Index(range(len(new_df['Well'].values)))],new_df) ##reindexing and keeping sample name as column
new_df = new_df[["Samples", "CT"]]
# print(new_df["Samples"].values)
# print(new_df)
# print(df_qubit.head(50))

#################################
## Sorting
#################################

index_numbers=[0]
arrayToBeFiltered=new_df['Samples'].values

for string in arrayToBeFiltered[1:]:
    m=re.search('(_)\w+', string) ## regex filter
    # print(m.group(0)[1:])
    index_numbers.append(int(m.group(0)[1:]))


new_df['index_numbers']=index_numbers
# print(new_df) ## checking if the index number column corresponds to the Samples column Tn5 number
sorted_new_df=new_df.sort_values(by=['index_numbers']).drop(['index_numbers'],axis=1)
# print(sorted_new_df)
sorted_new_df=sorted_new_df.set_index([pd.Index(range(len(sorted_new_df['Samples'].values)))],sorted_new_df)
# print(sorted_new_df) ## setting new index from 0 to len(Samples.values) -> 97
new_df = sorted_new_df
# print(new_df)

###############################
## Merge df_qubit and new_df dataframes and calculate correlation between CT and dsDNA conc
###############################

new_df = pd.merge(new_df, df_qubit, on=["Samples"]) ## Positive control lost on merging
mean_dna = new_df["dsDNA (ng/uL)"].mean()
# print(mean_dna)
correlation_CT_dna, p_val_CT_dna = scipy.stats.pearsonr(new_df["CT"], new_df["dsDNA (ng/uL)"])
# print("Correlation between CT value and dsDNA concentration:", correlation_CT_dna)
rsq_CT_dna = correlation_CT_dna ** 2
print("R squared value between CT values and dsDNA concentration:", rsq_CT_dna)
print("p value between CT values and dsDNA concentration:", p_val_CT_dna)

###################################
## Calculating adjusted_CT
####################################

# print(type(new_df["dsDNA (ng/uL)"]))
new_df["adjusted_CT"] = new_df["CT"] + np.log2(new_df["dsDNA (ng/uL)"]/mean_dna)
# print(new_df)
correlation_CTs, p_val_CTs = scipy.stats.pearsonr(new_df["CT"], new_df["adjusted_CT"])
rsq_CTs = correlation_CTs ** 2
print("R squared value bewtwen raw and adjusted CT values:", rsq_CTs)
print("p value between raw and adjusted CT values:", p_val_CTs)

correlation_adjCT_dna, p_val_adjCT_dna = scipy.stats.pearsonr(new_df["adjusted_CT"], new_df["dsDNA (ng/uL)"])
print("R squared value between adjusted_CT and dsDNA concentration:", correlation_adjCT_dna ** 2)
print("p value between adjusted_CT and dsDNA concentration:", p_val_adjCT_dna)
##############################
## Calculting the deltaCT
#############################

new_df["Positive"] = new_df.iloc[0]["adjusted_CT"]
new_df["deltaCT"] = new_df["adjusted_CT"] - new_df["Positive"]
new_df = new_df.drop(["Positive"], axis=1)
# print(new_df)
# print(new_df.head(50))
# print(new_df.loc["CT", 0])

##############################
## Calculating % digested
##############################

new_df["Percent_digested"] = (1 - (1/2**(new_df["deltaCT"]))) * 100
print("Mean pUC19 digested (calculated using adjusted_CT):", new_df["Percent_digested"].mean())
print(new_df)
new_df.to_csv(r"/Users/bholalamichhane/Desktop/Honours/Bioinformatics/tn5_qpcr_analysed_updated.csv",index=False)

#############################
## Plotting
##############################

# ## CT and adjusted_CT scatterplot
# x = new_df["CT"].values
# y = new_df["adjusted_CT"].values
# x = x.reshape(-1, 1)
# y = y.reshape(-1, 1)
# model2 = LinearRegression()
# model2.fit(x, y)
# plt.scatter(x, y)
# plt.plot(x, model2.predict(x), color="k")
# plt.xlabel("CT value")
# plt.ylabel("Adjusted CT value")
# plt.title("Correlation between CT and adjusted CT values")
# plt.text(13, 11, "rsq = 0.659")
# plt.text(13, 10.5, "p = 0.00")
#
# ## CT and dsDNA conc scatterplot
# x = new_df["CT"].values
# y = new_df["dsDNA (ng/uL)"].values # 10, not 9, so the fit isn't perfect
# x=x.reshape(-1, 1)
# y=y.reshape(-1, 1)
# model2 = LinearRegression()
# model2.fit(x, y)
# plt.scatter(x, y,color='g')
# plt.plot(x, model2.predict(x),color='k')
# plt.xlabel("Mean CT Value of each transposome")
# plt.ylabel("dsDNA concentration (ng/uL)")
# plt.title("Correlation between Ct Mean Value and DNA concentration")
# plt.text(13, 1.4, "rsq = 0.057")
# plt.text(13, 1.3, "p = 0.01")

# ## adjusted_CT and dsDNA conc scatterplot
# x = new_df["adjusted_CT"].values
# y = new_df["dsDNA (ng/uL)"].values # 10, not 9, so the fit isn't perfect
# x=x.reshape(-1, 1)
# y=y.reshape(-1, 1)
# model2 = LinearRegression()
# model2.fit(x, y)
# plt.scatter(x, y,color='g')
# plt.plot(x, model2.predict(x),color='k')
# plt.xlabel("Adjusted CT Value of each transposome")
# plt.ylabel("dsDNA concentration (ng/uL)")
# plt.title("Correlation between adjusted CT Mean Value and DNA concentration")
# plt.text(13, 1.4, "rsq = 0.113")
# plt.text(13, 1.3, "p = 0.001")
#
## bar graph of % pUC19 digested using adjusted_CT values
x_object = new_df["Samples"][1:]
# print(x_object)
y_pos = np.arange(len(x_object))
plt.bar(y_pos, new_df["Percent_digested"][1:], align='center', alpha=0.5, color="black")
plt.xticks(y_pos, new_df["Samples"][1:], rotation = "90")
plt.ylabel('% pUC19 digested')
plt.title('% pUC19 digested across Tn5 transposomes')
plt.tight_layout()
plt.style.use("ggplot")
plt.show()
