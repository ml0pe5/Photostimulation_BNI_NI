# Photostimulation_BNI_NI
Code related to the article "A computational biomarker of photosensitive epilepsy from interictal EEG"

The code provided in this repository can be used to compute functional networks from EEG data and to then place a computer model of ictogenicity (the theta model) on the networks. By using the theta model, the user can estimate the Brain Network Ictogenicity (BNI) and consequently the Node Ictogenicity (upon removal of nodes from a network). 

Two main functions are provided: the "plv_method" to compute the functional networks, and the "theta_model" to calculate the BNI. 
The "plv_method" uses two functions available in the "Auxiliary functions" directory.

All parameters and further details can be found in the article:
Lopes et al., "A computational biomarker of photosensitive epilepsy from interictal EEG"
