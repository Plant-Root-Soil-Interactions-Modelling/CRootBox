# CRootBox for rice soil column

This is the work from Trung-Hieu Mai, 2019 on nutrient uptake of rice root system. One part of the study was to reconstruct the rice root system based on root information from the rice soil column experiment. The root characteristics, such as the number of the nodal root, root radius, interbranch distances, were used as the input parameters. The total root mass and the root mass fractions between soil layers were used to fit the CRootBox model. For that reason, three new feature was implemented in the CRootBox: (1) non-linear function of the emergence of the nodal root, (2) distribution of different lateral types (S-type and L-type) along soil depth using a probability function, (3) scaling function of the interbranch distance between laterals along the soil depth.

To download and compile the work

```
git clone -b pubMai2019 https://github.com/Plant-Root-Soil-Interactions-Modelling/CRootBox.git

cd CRootBox

cmake.

make
```

then run the python scripts which generate root system and compare root mass fraction visually by ploting with 1:1 line 
```
python3 Mai2019_riceNERICA4_NoP_DryingV5.py
```

# More information about the study

see here [A model-data integration study for soil rice column using multiscale modelling approach considering rhizosphere gradients](https://www.researchgate.net/publication/331773573_A_model-data_integration_study_for_soil_rice_column_using_multiscale_modelling_approach_considering_rhizosphere_gradients)


