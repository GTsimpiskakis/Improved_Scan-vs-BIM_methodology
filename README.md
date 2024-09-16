# Face-vs-Segment comparison: a new method for comparing as-planned and as-built models

Research: Master's thesis

## Introduction
During my master’s thesis at TU/e, I developed an improved Scan-vs-BIM methodology, inspired by the work of **Jiang et al.** To test my idea, I implemented a Python algorithm. Although civil engineers are often reluctant to engage with coding, I embraced the challenge, and the results of my thesis research can be found here: [link]. Below, I outline the main steps of my research grouped into two big sections **RANSAC-DBSCAN segmentation method**, and **Face-vs-Segment Comparison method**, illustrating the two new methods proposed by this research.

## Table of Contents
**RANSAC-DBSCAN segmentation method**

1. Convert LAS to PLY
2. Point cloud denoising
3. Point cloud filtering, Downsampling, and RANSAC segmentation
4. Classification of RANSAC segments
5. Apply the DBSCAN algorithm on the wrong RANSAC segments
6. Apply the DBSCAN algorithm a second time to further split a segment
7. Remove noise from a correct DBSCAN Segment
8. Merge DBSCAN Segment
9. Store the names and the colours of each segment to a CSV file
10. Give unique colours to every Final segment
11. Merge all the segments in a single Point cloud model
12. Transformation procedure
13. Divide the transform point cloud model into Segments again.

**Face-vs-Segment Comparison method**

14. Parse BIM model in IFC file format
15. Adjust the area of faces in case of a Wall-wall face intersection
16. Comparison procedure between a Face and a Segment

## 1. Convert LAS to PLY: 
**Script name: 1_LAS to PLY_real.py**

The first step involves converting the raw point cloud model from the .las format to the .ply format using Script 1. This conversion is necessary to utilize the Open3D library in Python, which provides commands for point cloud processing of files in .ply format.

## 2. Point cloud denoising: 
**Script name: 2_SDK_real_Cube.py**

We apply the bounding box method to extract a desired part of the point cloud model. Since the point cloud model measured with laser scanner BLK360 represents an entire floor of a building under construction, we use the bounding box method to extract only the part of the floor that we need. Because the initial point cloud contains significant noise (Opposite buildings, very distant points, etc) we apply this method to extract only the useful points. Finally, the processed file depicted in Figure ... contains only the extracted portion (space of the measured indoor building environment of a building under construction) of the initial point cloud model.

## 3. Point cloud filtering, Downsampling, and RANSAC segmentation: 
**Script name: 3_Pre-procesing_Segmentation_real_data (RANSAC algorithm).py**

In this step, we use Statistical outlier removal (SOR) and Voxelisation methods to reduce the size of the point cloud model. In the same script, we apply the RANSAC algorithm to segment the point cloud model into segments. The results are depicted in Figure... As we can easily see, the results are not desirable as the algorithm has not identified all the building element surfaces in the point cloud model as separate segments. After several tests with different RANSAC parameters, and several tests using the DBSCAN algorithm it was realised that there was a problem in the segmentation procedure of the Point cloud model of Figure... that represents an indoor environment of a building under construction. 

To resolve this, there were three potential approaches:
1.  Develop a new segmentation algorithm from scratch.
2.  Conduct further literature research to find another suitable algorithm.
3.  Develop a semi-automated process by combining the RANSAC and DBSCAN algorithms.

Due to limited expertise in computer science, and limited time, the first two options were abandoned. Thus, the thesis research was focused on developing a new semi-automated method called the **RANSCA-DBSCAN segmentation method**.

The next chapters explain further the RANSCA-DBSCAN segmentation method.

## 4. Classification of RANSAC segments
**Script name: 4_Pre-procesing_real_data_RANSAC_segments_classification.py**

In Chapter 3 we applied the RANSAC algorithm and we produced 44 segments that are stored in the folder **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments**. From now on these 44 segments are named as **RANSAC segments**. In this step, we run the script 4. With this script, we iterate through each RANSAC segment and label them into three categories **Discarded**, **Wrong**, or **Correct**. For each label the segments are classified into three different folders in the folder path **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments**

Label **Discarded** means that these segments are useless and we delete them
Label **Wrong** means that the segment represents incorrectly a building element surface
Label **Correct** means that the segment represents correctly a building element surface

So from now on we only deal with the wrong RANSAC segments

## 5. Apply the DBSCAN algorithm on the wrong RANSAC segments
**Script name: 5_DBSCAN_Real_data_1st_DBSCAN.py** 

Script 5 is applied to each wrong segment generated from the RANSAC process of steps 3 and 4, located in the folder: **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\Wrong_Segments**. This script runs the DBSCAN algorithm. The user specifies the DBSCAN parameters required for each RANSAC segment produced in step 3 (For further explanations read the thesis document). The results are visualised in the Python environment and the user decides if the segmentation of the RANSAC segment is satisfactory.

the segmentation is unsatisfactory, the user adjusts the DBSCAN parameters and re-runs the script for the same RANSAC segment. Once the segmentation is deemed correct, the segments are imported into Meshlab or CloudCompare for further inspection. The user manually classifies the resulting DBSCAN sub-segments into two categories: **Correct** or **Discarded** (i.e., useless segments), which are stored in separate folders.

**Example:**
The Segment: **PLY_SDK_EU_Real_Segment_4_Wrong.ply** in the Wrong_Segments folder is illustrated in Figure...
![My Image Description](images/RANSAC_segment_4.png)

After applying script 5 with input DBSCAN parameters eps = 0.1 and MinPts = 5 Segment 4 is divided into six Segments as they are stored in folder path: **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\DBSCAN_clusters\PLY_SDK_EU_Real_Segment_4**. Then the user manually classifies these segments again into **Correct** and **Discarded**, to keep the correct segments in a separate file.

## 6. Apply the DBSCAN algorithm a second time to further split a segment

**Script name: 6_DBSCAN._SDK_Real_data_2nd_DBSCAN_segmentation (DBSCAN run for each seperate segment).py** 

In some cases, a segment generated from step 5 may require further segmentation. Script 6 is used to apply DBSCAN a second time to further divide these segments. For instance, Segment 24 located in 
**Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\DBSCAN_clusters\PLY_SDK_EU_Real_Segment_24** may need additional splitting.

## 7. Remove noise from a correct DBSCAN Segment

**Script name: 7_SDK_DBSCAN_Segment_Cube_Denoise.py**


After running the DBSCAN algorithm, the building element surfaces are represented by distinct segments. However, certain segments may still contain noise, which can distort the dimensions of the building elements. For instance, Segment 8_1 (located at **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\DBSCAN_clusters\PLY_SDK_EU_Real_Segment_8\Correct_Segments\PLY_SDK_EU_Real_Segment_8_1.ply**) has points on its left side that deviate from the actual edge of the segment, as shown in Figure...

To remove this noise, the bounding box method is applied. First, the segment is imported into CloudCompare, where the user manually identifies the coordinates of the bounding box that contains the noisy points. These coordinates are then provided as input to script 7, which deletes the unwanted points from the point cloud model. The final denoised segment (Segment 8_1_Denoised) is saved at: **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\DBSCAN_clusters\PLY_SDK_EU_Real_Segment_8\Correct_Segments\PLY_SDK_EU_Real_Segment_8_1_Denoised.ply**. The denoised version is depicted in Figure...

## 8. Merge DBSCAN Segment

**Script name: 8_Pre-procesing_Segments_Merging_ply_files_From DBSCAN_to_Real_segments.py**

After applying the DBSCAN algorithm to the RANSAC segments, it may be necessary to merge two or more DBSCAN segments to create a correct segment that accurately represents a building element surface. Script 8 is used to merge multiple point cloud segments (.ply files) together.

Example:...

## Final correct Segments
At this stage, we have identified all segments that best represent the building element surfaces in the indoor environment. These segments are stored in the following file path: **C:\Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\Final_Real_Segments\Attempt_1_Common_colours**.

## 9. Store the names and the colours of each segment to a CSV file

**Script name: 9_Pre-procesing_Segment_names_and_colours_in_CSV.py**

This script captures all the final segments produced by the RANSAC-DBSCAN segmentation process and stores their names and associated colors in a CSV file. This helps in identifying and managing each segment efficiently.

## 10. Give unique colours to every Final segment

**Script name: 10_Pre-procesing_Segments_Merging_ply_changing_Real_Segments colour.py**

Using the CSV file generated in the previous step, we manually identify segments that share the same color and store them in separate folders (e.g., **C:\Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\Final_Real_Segments\Change their colour**). Then by using Script 10 we give unique colours to each group of segments that share the same colour. The segments with the unique colours are stored in a different folder (ex. **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\Final_Real_Segments\Change their colour\Colour_148_103_189\Colour_148_103_189_changed**). Finally, each segment has a unique colour.

Finally we manually gather all the segments in the same folder **Python Output\Real Data\PLY files\SDK_real_dataset_1\Segments\Final_Real_Segments\Attempt_2_Unique_Colours**.

## 11. Merge all the segments in a single Point cloud model

**Script name: 11_Pre-procesing_Segments_Merging_Final_segments.py**

The segments from step 10 are merged into a single point cloud model. This consolidated point cloud will be used in the next step, which involves transforming the model to align with the BIM model.

## 12. Transformation procedure

**Script name: 12_GT_Topoplan_psich_SDK_real_data_Transformation(Translation_Rotation).py**

To compare the point cloud model with the BIM model, the point cloud needs to be aligned with the BIM model's coordinate system. Script 12 performs this crucial transformation step. Preferably, if there was a way to export the Point cloud model from the laser scanner device directly on the BIM model's Cartesian Coordinate System, the two models would be automatically aligned. However, due to the lack of this capability, e manually overlap the models by transforming the point cloud to match the BIM location. 

The transformation process is executed in two steps:
1.  **Translation**: Positioning the point cloud in the correct location.
2.  **Rotation**: Adjusting the orientation to match the BIM model. (For more details, refer to the thesis document.)


## 13. Divide the transform point cloud model into Segments again.

**Script name: 13_Pre-procesing_Segments_Devide_the mergedPCM_into_Segments_again.py**

After we have transformed the Point cloud model into the location of the BIM model, the point coordinates of the point cloud have been changed. The next step is to divide the point cloud model into the same segments asas those produced in step 10. Since each segment is uniquely coloured, we can use this feature to recreate the segments as they were. The final segments after this step are saved in the folder: **Python Output\Real Data\PLY files\SDK_real_dataset_final\Segments**.


## Face-vs-Segment Comparison method

**Script name: 14 Face-vs-Segment_comparison_method.py**

This research introduces an improved Scan-vs-BIM methodology by proposing a novel comparison method between a point cloud model and a BIM model, called the **Face-vs-Segment comparison method**. This new approach compares building elements at the surface level, offering a more detailed and precise way of comparing the two models. This new approach compares building elements at the surface level, offering a more detailed and precise way of comparing the two models. 

The key innovation of this method is that it compares a **Face** (a building element surface in the BIM model) with a **Segment** (the corresponding surface in the point cloud model). By doing so, the method ensures a more accurate and detailed comparison between the two models.

Below is an overview of the key steps involved in developing the Face-vs-Segment comparison method, referencing specific lines of code from script 14 responsible for each step.

## 14. Parse BIM model in IFC file format
The first step in the Face-vs-Segment comparison method is to parse the IFC file of the BIM model through its classes and entities to identify the coordinates of all the building elements that are included in this analysis (such as walls, columns, slabs, beams, and opening elements). The different orientations of the the building elements in the BIM model as well as the IFC type of building elements can indicate different types of geometric representations in the IFC schema (See thesis document for more info). 

Thus, lines 12-1163 of script 14 are dedicated to extracting and identifying the coordinates of each surface (face) of the building elements.

Each building element (e.g., a wall or slab) has six **Faces**. By the end of line 1163, the script generates a **complete matrix** for each building element, where each row contains six arrays representing the four corner coordinates of each face. This matrix also stores, in each row, metadata about the building elements, such as their IFC types and geometric details.

## 15. Adjust the area of faces in case of a Wall-wall face intersection

Lines 1248-1442 explain the case in which we have a wall-wall intersection. In this case, how the BIM model has been modelled by the modeller matters. Figure... illustrates a floor plan of a wall-wall intersection. In the case of a wall-wall intersection, the wall face of one wall can be entirely folded inside the wall face of the other wall. This can lead to inaccurate results in the later steps of comparing the wall face with the corresponding point cloud segment. For this reason, we deduct the area of the face that is entirely folded inside the second wall face, from the are of the second wall face.

## 16. Comparison procedure between a Face and a Segment

The rest of the code deals with comparing a face with a segment.  Several conditions and Rules are formulated to govern this comparison procedure. As it is difficult to explain thoroughly this process in a readme file you can also search in the thesis Document for further info. 

## Conclusion

The Face-vs-Segment comparison method offers a more refined approach to comparing point cloud and BIM models, focusing on building element surfaces. This method introduces a detailed comparison at the face and segment level, enhancing accuracy in Scan-vs-BIM applications.

I hope this research provides valuable insights and improves current methodologies! Feel free to contact me if you have any questions.




