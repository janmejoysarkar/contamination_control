
![Logo](https://suit.iucaa.in/sites/default/files/top_banner_compressed_2_1.png)


# SUIT Contamination Correction 🧹 ☀️ 🛰️ 

The SUIT instrument on Aditya-L1 maintains its image sensor at a very low temperature (-55 degC). This causes some volatiles to condense on the CCD surface. Traces of these contaminants are seen in the images recorded by SUIT.
This collection of modules is necessary to remove the contaminants from the SUIT images.


- `full_disk_flat_gen.py`: Generates calibration file. Saves in `data/interim`
- `full_disk_process.py`: Processes files in `data/raw` using corresponding calibration file in `data/interim`
- `roi_correction.py`: To apply contamination correction for RoI Images. Works best for NB03 and NB04 feature rich images.
- `validation_with_iris.py`: Used to validate photometry with IRIS SJI images.

## Screenshots

![Correction of contaminants on SUIT NB05 filter](README_files/Figure_1.png)

## Authors

- [@janmejoysarkar](https://github.com/janmejoysarkar)

## Acknowledgements

 - [ISRO, Aditya-L1](https://www.isro.gov.in/Aditya_L1.html)
