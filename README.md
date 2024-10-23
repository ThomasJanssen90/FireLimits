# FireLimits
This repository contains the code for reproducing the results from the paper "What are the limits to the growth of boreal fires?" by Thomas Janssen and Sander Veraverbeke. The code for the main analysis is found [here](Analysis/Analysis_firebarriers.R). The specific analysis scripts that are dependencies for the main code are found [here](Analysis/SourceFunctions).

## Abstract
Several boreal forest regions, including East Siberia, experienced elevated fire activity in recent years. These fires substantially contributed to greenhouse gas emissions and air pollution. We currently have incomplete understanding of what eventually stops the spread and thus limits the growth of these fires. This hampers our ability to model their dynamics and predict their impacts. Here, we investigated the locations and timing of fire stops across the vast taiga of East Siberia.We used several geospatial datasets, including fire weather and landscape variables, which were overlaid with fire perimeters recorded between 2012 and 2022 to characterize the causes of fire stops. Using this approach, we attributed 87% of all fire stops to a statistically significant (p < 0.01) change in one or more of these environmental variables over time (fire weather) or space (landscape).

## Study area
Our analysis considered the fire prone eastern Siberian boreal forest, also referred to as taiga, and mountain tundra ecoregions. The extent of our study domain was delineated based on ecoregions, geographic areas with similar ecological and environmental characteristics (Olson et al., 2001). We selected the East- and Northeast Siberian taiga ecoregions as well as the adjacent Cherskii-Kolyma Mountain tundra and the Trans-Baikal Bald Mountain tundra ecoregions (Fig. 1). 

![Figure1-01](https://github.com/user-attachments/assets/e0e50911-9d51-4828-a447-a8dd6761ed2a)
Figure 1 a) ecoregions, b) land cover, c) permafrost landscapes and d) fire ignition sources in the study domain. The grey lines denote country boundaries and the black lines ecoregion boundaries. Ecoregions are derived from Olson et al. (2001), land cover is adapted from Bartalev et al. (2003, 2010), permafrost distribution is derived from Obu et al. (2019), Yedoma presence is from Strauss et al. (2022) and ignition sources are from Janssen et al. (2023). The landcover and permafrost distribution were both aggregated by majority class to a 5 by 5 km spatial resolution for easier visual interpretation, the ignition source data was maintained at its original 0.5-degree spatial resolution. See Methods for a detailed description of the datasets.

## Obtaining fire perimeter locations.
The focus of our analyses is the fire perimeter location, a single spatial point along the fire perimeter to which we can assign a potential cause of the fire stopping there at a specific time (Fig. 2). Before we obtained these locations, we omitted fires smaller than six contiguous pixels (0.54 km², which was also the median fire size) from the dataset because our focus was on large wildfires. Subsequently, for each year in the harmonized burned area dataset (2012-2022), we buffered the burned area pixels (300 m) by a single pixel following first order Queen’s adjacency so that a -at least- 300-meter-wide buffer surrounded all observed fire pixels. These pixels were then converted to spatial point locations with an x and y coordinate, the fire perimeter locations (Fig. 2). We then extracted the burn date of the nearest burned pixel to obtain the last known burn date at the perimeter. The last know burn date is assumed to be the same or very close to the day that the fire stopped spreading at that location. For efficient data management and organization, each perimeter location was assigned a unique identifier. This procedure yielded on average around 227.000 perimeter locations per year over the entire study area, with up to 346.000 perimeter locations in the record year 2020.

![Figure2](https://github.com/user-attachments/assets/c5009033-5775-432e-abd2-38ea246d6a14)
Figure 2 Attributing top-down climate and bottom-up landscape drivers impeding fire growth in a subset of the study domain. a) timeseries of vapor pressure deficit (VPD) surrounding the last recorded burn date at a single perimeter location showing a structural break (dashed line) and a decline in VPD at the last recorded burn date. The grey shaded area depicts the one day before and two days after the last recorded burn date period in which a sudden decline in VPD can result in VPD being flagged as potential cause. b) the harmonized burned area data of the year 2020 and the perimeter locations. c) establishing the direction perpendicular to the fire perimeter by identifying the nearest burned pixels d) using the direction perpendicular to the fire perimeter to determine the orientation of the ellipse e) creating a 100x100 meter regular grid of burned and unburned locations inside the polygon to extract landscape variables f) extracting the global surface water values to the grid. g) the global surface water data and the perimeter locations with surface water attribution.

## Attributing limitations to fire spread.
The perimeter locations and their associated last recorded burn date indicate the location and timing of fire stops, yet not their potential causes. To attribute the causes of the fire stops we carried out two types of analyses. In these two analyses, every separate potential driver could be flagged as a cause so that multiple causes for a single perimeter point were possible.

The first analysis was carried out by using the high temporal resolution of the ERA5-Land derived meteorological variables in a time series analysis. At each perimeter point, we extracted the hourly vapor pressure deficit (VPD), [code link VPD](Analysis/SourceFunctions/VPD_Functions.R), wind fire spread index (WFSI), [code link WFSI](Analysis/SourceFunctions/Wind_Functions.R), and surface soil moisture data, [code link soil moisture](Analysis/SourceFunctions/SoilMoisture_Functions.R), for the 10 days surrounding the last recorded burn date (Fig. 2a). For all three variables we first performed a Students t-Test using the R package stats to establish whether there was a significant (p < 0.01) difference between the five days before and including the last recorded burn date and the five days after the last recorded burn date. If a significant decline in VPD, WFSI or a significant increase in soil moisture was observed that could indicate a cause for the fire stop, we proceeded to the next step. In the next step, to detect the presence of a significant sudden change or break point in the time series, we used the method developed by Bai and Perron (2003) and implemented in the R package strucchange to derive single or multiple break points in the time series. If a breakpoint was detected, and if the breakpoint (including its confidence interval) was timed one day before to two days after the last recorded burn date, we flagged the perimeter point as potentially caused by a change in the specific weather variable (Fig. 2a).

The second analysis leveraged the high spatial resolution of the landscape variables, [code link landscape drivers](Analysis/SourceFunctions/Landscape_Functions.R), to establish whether there was a significant change of a variable in space along a gradient from inside to outside of the fire perimeter. For this analysis, at each fire perimeter point, an elliptical polygon was created oriented with its major axis perpendicular to the fire perimeter (Fig. 2b-2d). The major axis length was set to three times and the minor axis length to one and a half time the spacing between the fire perimeter points, which is the same as the spatial resolution of the harmonized burned area data (300 m). The centre of the ellipse was set halfway the perimeter point location (the centroid of the first non-burned pixel) and the nearest burned pixels, so that approximately half of the area of the ellipse overlapped with burned pixels and the other half overlapped with unburned pixels (Fig. 2e). 

Subsequently, the data for all landscape variables were extracted to the ellipse at a 100 m resolution spacing, a spatial resolution in between the highest resolution landscape datasets (30 m) and the lowest spatial resolution land cover map (345 m, Table 1). In that way, we obtained two groups with observations of landscape variables representing to halves of the elliptical polygon, a group of burned locations inside the perimeter and a group of unburned locations located outside the perimeter (Fig. 2f). We performed two types of tests to establish a statistically significant difference in the landscape characteristics between these groups. For continuous variables such as above-ground biomass (AGB), percentage tree cover and burn history, we performed a t-Test. For the categorical variables such as land cover type and the presence of roads and water, we performed a Χ2 test. In case of a statistically significant difference (p < 0.01) between the groups for road and water presence, we also confirmed that the road and water presence were in fact higher outside the perimeter by comparing the means of the groups. Similarly, we confirmed that if a significant land cover change was detected, that the relatively highest increase in the area cover of a land cover class outside the perimeter was different than the most common land cover class inside the fire perimeter. Finally, for the downslope we tested whether there was a significant decline in elevation from inside to outside the perimeter and only flagged perimeter locations as being caused by downslope if the slope was lower or equal to -5°, which is a threshold value found in experimental studies below which downslope substantially limits fire spread (Sullivan et al., 2014). 

In our analysis, multiple potential causes can be attributed to a single perimeter location. It is also expected that multiple drivers often collectively contribute to fire quenching, for example when fire barriers such as roads and rivers are only effective if weather conditions change simultaneously (e.g. Holsinger et al., 2016; Moreno et al., 2014). However, our list of potential drivers included some drivers that supersede the effect of other drivers. For example, if a fire is extinguished at the banks of a large river, in our analysis this location will be flagged as caused by a change in surface water, a change in fuel load and a change in land cover. To prevent this double counting, we used a simple hierarchical decision tree to separate these interdependent drivers (Fig. 3). The two fire weather drivers; an increase in fuel moisture and a decline in the WFSI, could always be flagged separately as they can be considered independent from one and other and independent of the other drivers. Also, the presence of a pronounced downslope in the landscape is independently flagged as a potential cause (Fig. 3). The remaining potential fire stop drivers were all related to the availability of fuels, and therefore these were hierarchically ordered going from relatively hard fire barriers such as water bodies to softer landscape changes. The presence of surface water (1) and roads (2) were considered hard barriers, offering no available fuel. Hereafter, the burn history (3) was chosen in the decision tree as a previous fire would have substantially impacted the availability of fuels (4) as well as the land cover (5). The availability of fuels (4) was chosen over a change in land cover (5) as a decline in fuel load is directional while for a change in land cover, we only test for a significant change in land cover and not a change to a land cover less conductive to fire spread per se. In this way, at each perimeter location a maximum total of four independent drivers could be flagged as a potential cause of impeding fire spread at the perimeter. Because of this, for most statistics and figures we used the percentage of perimeter points where an independent driver was flagged which therefore could add up to a total of more than 100%. For the figures showing relative contributions (e.g. stack bars) where the total must add up to 100%, we divided the contribution of every driver by the total number of flagged drivers at each perimeter point.

![Methods_flowdiagram-01](https://github.com/user-attachments/assets/886860b3-efba-44b3-934c-723058dca72a)
Figure 3 Attributing changing weather and landscape drivers as potential limits to fire growth using a decision tree model. The effect of a change in the climate drivers related to fuel moisture and the wind fire spread index (WFSI) were tested in the form of a time series analysis. The effect of a change in landscape drivers were tested separately as a change in these variables over space from inside to outside the fire perimeter (see Fig. 2).  

## Proof of concept
To provide a proof of concept of our method, a small subset of the study area was selected for a detailed study on how the algorithm performed on two average sized fires that burned for about one month in 2021 (Fig. 4). Because of limited cloud cover in this period, we were able to use six high quality Sentinel-2 images to track and visualize the fire spread and the development of the fire perimeter. At the end of July 2021, the most southern fire is spreading while confined by the Lena River in the south and a smaller tributary in the Northwest (Fig. 4a-c). On the 8th of August, the outer perimeter of the southernmost fire is almost complete, while the northern fire is still spreading (Fig. 4d). Finally, on the 25th of August, both fires have been completely quenched and no new perimeter is added (Fig. 4f). While the distance between the outer perimeters of both fires is only 5 km, the attributed drivers impeding fire spread were very different. In the early southern fire, most of the perimeter was attributed to bottom-up landscape drivers, mostly the presence of surface water, and to some degree a decline in fuel load and land cover change (Fig. 4f). This in contrast to the later northern fire that was largely quenched by a top-down increase in fuel moisture, that is a decline in VPD or increase in soil moisture (Fig. 4f). This transition can also be observed in the time series of fire perimeter attribution (Fig. 4g), which starts off as heavily dominated by the presence of surface water before the 2nd of August, followed by a mix of landscape drivers and a decline in the WFSI and after the 8th of August by fuel moisture as VPD declines below the threshold value of 1.2 kPa (Fig. 4e). In general, despite the difference in spatial resolution between the harmonized burned area data (300 m) from which the perimeters were derived and the Sentinel-2 images (20 m), the algorithm performs well in ascribing large parts of the fire perimeter to the presence of water bodies, both in the form of rivers and streams in the southern fire as well as solitary lakes in the northern fire. 

![Figure4Sentinel](https://github.com/user-attachments/assets/65886ae4-89ae-4d71-892d-b08f549ec715)
Figure 4 A demonstration of the methodology on two averaged size fires that burned in 2021 along the Lena River, 350 km northwest of Yakutsk, Russia. Panels a to f show the progression of burning using a Seninel-2 false colour composite to visualize the active fires, smoke plumes and burned area. The acquisition dates are in the top of the panel, as well as an arrow representing the wind direction and speed derived from ERA5-Land hourly reanalysis data. In panels a to f, the perimeter locations on and before the date of image acquisition are plotted with the colours representing to which subset of potential drivers the point was attributed. g) the daily sum of the perimeter attributed to a subset of potential drivers and e) the time series of average vapor pressure deficit at 2 m in the area derived from ERA5-Land reanalysis data with the image acquisition dates indicate with the downward pointing arrows and associated panel letters. The grey horizontal arrow indicates the 1.2 kPa VPD threshold below which fire activity is found to be extremely rare in these boreal ecosystems, derived from Balch et al. (2022) and Clarke et al. (2022). 

