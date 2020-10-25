#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
#Import file
mouse_metadata_path="C:/Users/sledg/PythonData/matplolib-challenge/Resources/mouse_metadata.csv"
study_results_path="C:/Users/sledg/PythonData/matplolib-challenge/Resources/study_results.csv"
#Read Results"
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)
#Merge Data
merge_data= pd.merge(mouse_metadata, study_results)
merge_data


# In[2]:


#Check number of mice
merge_data["Mouse ID"].count()


# In[3]:


#Deduplicate
clean_data=merge_data.drop_duplicates(subset=["Mouse ID", "Timepoint"], keep="last")
clean_data


# In[4]:


#Clean
clean_data["Mouse ID"].count()


# In[5]:


"Summary Statistics"
sum_stat_group = clean_data.groupby("Drug Regimen")
regimen_mean = pd.DataFrame(sum_stat_group["Tumor Volume (mm3)"].mean())
regimen_median = pd.DataFrame(sum_stat_group["Tumor Volume (mm3)"].median())
regimen_variance = pd.DataFrame(sum_stat_group["Tumor Volume (mm3)"].var())
regimen_stdev = pd.DataFrame(sum_stat_group["Tumor Volume (mm3)"].std())
regimen_SEM = pd.DataFrame(sum_stat_group["Tumor Volume (mm3)"].sem())
whole=regimen_mean.merge(regimen_median, how="outer", on="Drug Regimen", suffixes=("_mean", "_median"))
whole=whole.merge(regimen_variance, how="outer", on="Drug Regimen")
whole=whole.merge(regimen_stdev, how="outer", on="Drug Regimen")
whole=whole.merge(regimen_SEM, how="outer", on="Drug Regimen")
whole.rename(columns={"Tumor Volume (mm3)_x": "Variance", "Tumor Volume (mm3)_y": "stdev", "Tumor Volume (mm3)": "SEM"})


# In[6]:


"Bar and Pie Charts"
# bar chart pandas 
mice_treat = clean_data.groupby(["Drug Regimen"]).count()["Mouse ID"]
mice_treat.plot(kind="bar")
plt.xlabel("Drug Regimen",fontsize = 14)
plt.ylabel("Number Mice", fontsize= 14)
plt.legend(loc="best")
plt.title("Mice per Treatment", fontsize=14)
plt.show()


# In[7]:


#bar chart pyplot
mice_list =(clean_data.groupby(["Drug Regimen"])["Mouse ID"].count()).tolist()
mice_list

x_axis = np.arange(len(mice_treat))
plt.bar(x_axis, mice_list, color='b', alpha=0.8, align='center')
tick_locations = [value for value in x_axis]
plt.xticks(tick_locations, ['Capomulin', 'Ceftamin', 'Infubinol', 'Ketapril', 'Naftisol', 'Placebo', 'Propriva', 'Ramicane', 'Stelasyn', 'Zoniferol'],  rotation='vertical')
plt.xlim(-0.75, len(x_axis)-0.25)
plt.ylim(0, max(mice_list)+10)

plt.title("Mice Per Treatment",fontsize = 14)
plt.xlabel("Drug Regimen",fontsize = 14)
plt.ylabel("Number Mice",fontsize = 14)


# In[8]:


#pie chart gender df
groupby_gender = clean_data.groupby(["Mouse ID","Sex"])
groupby_gender
gender_df = pd.DataFrame(groupby_gender.size())
mouse_gender_df = pd.DataFrame(gender_df.groupby(["Sex"]).count())
mouse_gender_df.columns = ["Total Count"]

# create and format the percentage of female vs male
mouse_gender_df["Percentage of Sex"] = (100*(mouse_gender_df["Total Count"]/mouse_gender_df["Total Count"].sum()))

# format the "Percentage of Sex" column
mouse_gender_df["Percentage of Sex"] = mouse_gender_df["Percentage of Sex"]

# gender_df
mouse_gender_df


# In[9]:


#pie chart pandas
colors = ["green", "blue"]
explode = (0.1, 0)
plot = mouse_gender_df.plot.pie(y="Total Count",figsize=(12,10), colors = colors, startangle=140, explode = explode, shadow = True, autopct="%1.1f%%")
plt.title("Gender Mouse Population",fontsize = 14)
plt.ylabel("Sex",fontsize = 14)
plt.axis("Equal",fontsize = 14)


# In[10]:


#pie chart pyplot
labels = ["Female","Male"]
sizes = [49.8,50.2]
colors = ["g", "b"]
explode = (0.1, 0)
plt.pie(sizes, explode=explode,labels=labels, colors=colors, autopct="%1.1f%%", shadow=True, startangle=140,)
plt.title("Gender Mouse Population",fontsize = 14)
plt.ylabel("Sex",fontsize = 14)
plt.axis("Equal",fontsize = 14)


# In[11]:


"Quartiles, Outliers, and Boxplots"
# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin List
Capomulin_df = clean_data.loc[clean_data["Drug Regimen"] == "Capomulin",:]
Capomulin_last = Capomulin_df.groupby('Mouse ID').max()['Timepoint']
Capomulin_vol = pd.DataFrame(Capomulin_last)
Capomulin_merge = pd.merge(Capomulin_vol, clean_data, on=("Mouse ID","Timepoint"),how="left")
Capomulin_merge.head()


# In[12]:


#Capomulin Print
Capomulin_tumors = Capomulin_merge["Tumor Volume (mm3)"]
quartiles =Capomulin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq
lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Lower quartile of Capomulin tumors= {lowerq}")
print(f"Upper quartile of Capomulin tumors= {upperq}")
print(f"Interquartile range of Capomulin tumors= {iqr}")
print(f"Median of Capomulin tumors= {quartiles[0.5]} ")
print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[13]:


#Ramicane List
Ramicane_df = clean_data.loc[clean_data["Drug Regimen"] == "Ramicane", :]
Ramicane_last = Ramicane_df.groupby('Mouse ID').max()['Timepoint']
Ramicane_vol = pd.DataFrame(Ramicane_last)
Ramicane_merge = pd.merge(Ramicane_vol, clean_data, on=("Mouse ID","Timepoint"),how="left")
Ramicane_merge.head()


# In[14]:


# Ramicane Print
Ramicane_tumors = Ramicane_merge["Tumor Volume (mm3)"]
quartiles =Ramicane_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq
lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Lower quartile of Ramicane tumors= {lowerq}")
print(f"Upper quartile of Ramicane tumors= {upperq}")
print(f"Interquartile range of Ramicane tumors= {iqr}")
print(f"Median of Ramicane tumors= {quartiles[0.5]} ")
print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[15]:


#Infubinol List
Infubinol_df = clean_data.loc[clean_data["Drug Regimen"] == "Infubinol", :]
Infubinol_last = Infubinol_df.groupby('Mouse ID').max()['Timepoint']
Infubinol_vol = pd.DataFrame(Infubinol_last)
Infubinol_merge = pd.merge(Infubinol_vol, clean_data, on=("Mouse ID","Timepoint"),how="left")
Infubinol_merge.head()


# In[16]:


#Infubinol Print
Infubinol_tumors = Infubinol_merge["Tumor Volume (mm3)"]
quartiles =Infubinol_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq
lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Lower quartile of Infubinol tumors= {lowerq}")
print(f"Upper quartile of Infubinol tumors= {upperq}")
print(f"Interquartile range of Infubinol tumors= {iqr}")
print(f"Median of Infubinol tumors= {quartiles[0.5]} ")
print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[17]:


#Ceftamin List
Ceftamin_df = clean_data.loc[clean_data["Drug Regimen"] == "Ceftamin", :]
Ceftamin_last = Ceftamin_df.groupby('Mouse ID').max()['Timepoint']
Ceftamin_vol = pd.DataFrame(Ceftamin_last)
Ceftamin_merge = pd.merge(Ceftamin_vol, clean_data, on=("Mouse ID","Timepoint"),how="left")
Ceftamin_merge.head()


# In[18]:


#Ceftamin Print
Ceftamin_tumors = Ceftamin_merge["Tumor Volume (mm3)"]
quartiles =Ceftamin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq
lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)
print(f"Lower quartile of Ceftamin tumors= {lowerq}")
print(f"Upper quartile of Ceftamin tumors= {upperq}")
print(f"Interquartile range of Ceftamin tumors= {iqr}")
print(f"Median of Ceftamin tumors= {quartiles[0.5]} ")
print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[19]:


#box plot
plot = [Capomulin_tumors, Ramicane_tumors, Infubinol_tumors, Ceftamin_tumors]
Regimen= ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']
fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.set_title('Tumor Volume For Each Mouse',fontsize =14)
ax1.set_ylabel('Final Tumor Volume (mm3)',fontsize = 14)
ax1.set_xlabel('Drug Regimen',fontsize = 14)
ax1.boxplot(plot, labels=Regimen, widths = 0.4, patch_artist=True,vert=True)


# In[37]:


"Line and Scatter Plot"
# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
Campomulin_line_df = Capomulin_df.loc[Capomulin_df["Mouse ID"] == "b128",:]
x_axis = Campomulin_line_df["Timepoint"]
tumorsize = Campomulin_line_df["Tumor Volume (mm3)"]
fig1, ax1 = plt.subplots(figsize=(10, 10))
plt.title('Capomulin treatmeant of mouse b128',fontsize =25)
plt.plot(x_axis, tumorsize,linewidth=3, markersize=15,marker="o",color="b")
plt.xlabel('Timepoint (Days)',fontsize =14)
plt.ylabel('Tumor Volume (mm3)',fontsize =14)
plt.show()


# In[21]:


## Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
fig1, ax1 = plt.subplots(figsize=(10, 10))
avg_tumor_vol=Capomulin_df.groupby(['Mouse ID']).mean()
marker_size=10
plt.scatter(avg_tumor_vol['Weight (g)'],avg_tumor_vol['Tumor Volume (mm3)'],s=175, color="b")
plt.title('Average Tumor Volume vs Mouse Weight',fontsize =14)
plt.xlabel('Weight (g)',fontsize =14)
plt.ylabel('Averag Tumor Volume (mm3)',fontsize =14)


# In[35]:


"Correlation and Regression"
# Calculate the correlation coefficient and linear regression model for mouse weight and average tumor volume for the Capomulin regimen
corr=round(stats.pearsonr(avg_tumor_vol['Weight (g)'],avg_tumor_vol['Tumor Volume (mm3)'])[0],2)
x_values = avg_tumor_vol['Weight (g)']
y_values = avg_tumor_vol['Tumor Volume (mm3)']
(slope, intercept, r_value, p_value, std_err) = stats.linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_equation = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))

fig1, ax1 = plt.subplots(figsize=(10, 10))
plt.scatter(x_values,y_values,s=175, color="b")
plt.plot(x_values,regress_values,"r-")
plt.title('Regression Plot Avg Tumor Volume vs Mouse Weight',fontsize =14)
plt.xlabel('Weight(g)',fontsize =14)
plt.ylabel('Average Tumore Volume (mm3)',fontsize =14)
ax1.annotate(line_equation, xy=(10, 20), xycoords='data', xytext= (0.5, 0.95),textcoords='axes fraction',horizontalalignment='right', verticalalignment='top',fontsize=14,color="g")

print(f"The r-squared is: {r_value**2}")


# In[ ]:




