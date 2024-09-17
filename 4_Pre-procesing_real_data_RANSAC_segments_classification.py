#B.2. Point cloud file Downsampling
#First this script executes ply file downsampling without any user's internvention

#B.3. Point cloud file Segmentation
#Second executes (denoised) file segmentation by applying RANSAC method
import numpy as np
import open3d as o3d
import matplotlib.pyplot as plt 
import os 
import subprocess
#pcd_file_path500k = r"C:/Users/dell/OneDrive/Υπολογιστής/CME/Season 2023-2024/Graduation Project/Thesis/Data/Test Data/Topoplan/Produced PCD from test code 1/Apartment.ply"
ply_file_path_in=  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segments_2.ply"
ply_file_path_out =  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_DownSampled.ply"
#pcd_file_path500k = r"C:/PCDenwseis1.pcd"
pcd_file_path_out =  r"C:/Python Output/Real Data/PCD files/SDK_real_dataset_1/PLY_SDK_EU_Real_DownSampled.pcd"

def read_pcd(pcd_file_path):
    pcdfile = o3d.io.read_point_cloud(pcd_file_path)
    return pcdfile
pcd = read_pcd(ply_file_path_in)
print(pcd)
 #Data Pre-processing
#pcd_center = pcd.get_center()
#print(pcd_center)
#pcd_Translated = pcd.translate(pcd_center)
#print(pcd_Translated)

#Statistical outlier filter
#nn = 8  #Number of nearest neighbours
#std_multiplier = 3
#filtered_pcd = pcd.remove_statistical_outlier(nn, std_multiplier)

#outliers = pcd.select_by_index(filtered_pcd[1], invert = True)
#outliers.paint_uniform_color([1,0,0])
#filtered_pcd = filtered_pcd[0]
#print(outliers)

#o3d.visualization.draw_geometries([filtered_pcd, outliers])
#Voxel_size= 0.01
#pcd_downsample = pcd.voxel_down_sample(voxel_size = Voxel_size)
#o3d.visualization.draw_geometries([pcd_downsample])
# Save Open3D PointCloud object to PCD file
#o3d.io.write_point_cloud(pcd_file_path_out, pcd_downsample, write_ascii=True)
#print(pcd_downsample)

def read_pcd2(pcd_file_path):
    pcd = o3d.io.read_point_cloud(pcd_file_path)
    points = np.asarray(pcd.points)
    colors = np.asarray(pcd.colors)
    return points, colors

    #create_TXT(txt_file_path,Sample_with_ID), Call command
def create_TXT(txt_Segement_path, segments):
    
    with open(txt_Segement_path, "w") as f:
        for i in range(len(segments)):
            segment_points = np.asarray(segments[i].points)
            segment_colors = np.asarray(segments[i].colors)
            for point, color in zip(segment_points, segment_colors): #iterate through point and color simultaneously in order to write each point with its color in the new txt file, and arrays
                f.write(f"{point[0]:.6f} {point[1]:.6f} {point[2]:.6f} {color[0]:.6f} {color[1]:.6f} {color[2]:.6f}\n")

def read_txt(txt_file_path):
    points = []
    colors = []
    with open(txt_file_path, "r") as f:
        for line in f:
            data = line.strip().split()
            point = [float(data[i]) for i in range(3)]
            color = [float(data[i]) for i in range(3, 6)]
            points.append(point)
            colors.append(color)
    return np.array(points), np.array(colors)
                                                      
def write_ply_from_pcd(points, colors, filename):
    vertex = np.zeros(points.shape[0], dtype=[("x", "f4"), ("y", "f4"), ("z", "f4"), ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    vertex["x"] = points[:, 0]
    vertex["y"] = points[:, 1]
    vertex["z"] = points[:, 2]
    vertex["red"] = colors[:, 0]
    vertex["green"] = colors[:, 1]
    vertex["blue"] = colors[:, 2]

    with open(filename, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        for i in range(len(points)):
            f.write("{:.3f} {:.3f} {:.3f} {:.0f} {:.0f} {:.0f}\n".format(points[i][0], points[i][1], points[i][2], colors[i][0], colors[i][1], colors[i][2]))
def write_ply_from_np_array(array, filename):
    vertex = np.zeros(points.shape[0], dtype=[("x", "f4"), ("y", "f4"), ("z", "f4"), ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    vertex["x"] = array[:, 0]
    vertex["y"] = array[:, 1]
    vertex["z"] = array[:, 2]
    vertex["red"] = array[:,3]*255
    vertex["green"] = array[:, 4]*255
    vertex["blue"] = array[:, 5]*255

    with open(filename, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        for i in range(len(points)):
            f.write("{:.3f} {:.3f} {:.3f} {:.0f} {:.0f} {:.0f}\n".format(array[i][0], array[i][1], array[i][2], array[i][3], array[i][4], array[i][5]))

points, colors= read_pcd2(pcd_file_path_out)
write_ply_from_pcd(points, colors *255, ply_file_path_out) #Scale colors back to the 0-255 range for PLY
print("Conversion completed.")
j = 0
if j == 1:
    #Multi-order RanSAC planar Segmentation
    max_plane_index =44
    plane_to_pnt_dist = 0.02
    segment_models = {}
    segments = {}
    rest = pcd
    Segments_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1//Segments/PLY_SDK_EU_Real.ply"
    for i in range(max_plane_index):
        colors = plt.get_cmap( "tab20")(i)
        segment_models[i], inliers = rest.segment_plane(distance_threshold = plane_to_pnt_dist, ransac_n = 3, num_iterations = 7000) 
        #we exclude inliers as indexes that can be interpreted only by open3d library's functions
        segments[i]= rest.select_by_index(inliers)
        segments[i].paint_uniform_color(list(colors[:3]))
        rest= rest.select_by_index(inliers, invert = True)
        directory, file_name_with_extension = os.path.split(Segments_path)
        file_name, extension = os.path.splitext(file_name_with_extension)
        additional_word = f"_Segment_{i+1}"
        # Modify the filename based on the condition
        new_file_name = file_name + additional_word
        # Recreate the output file paths
        Desirable_Segm_path_out = os.path.join(directory, new_file_name + extension)
        segment_points = np.asarray(segments[i].points)
        segment_colors = np.asarray(segments[i].colors)
        write_ply_from_pcd(segment_points, segment_colors*255, Desirable_Segm_path_out)

        print("pass", i,"/", max_plane_index, "done")
        #o3d.visualization.draw_geometries([inliers] )
    o3d.visualization.draw_geometries([segments[i]for i in range(max_plane_index)]+ [rest])
    #print(segments)
    #print(segment_models)
    txt_Segments_path = r"C:/Python Output/Real Data/TXT files//PLY_SDK_EU_Real_segements.txt"
    create_TXT(txt_Segments_path, segments)
    segm_points,segm_colors= read_txt(txt_Segments_path)
    Desireble_Segm_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segements_2.ply"
    write_ply_from_pcd(segm_points, segm_colors*255, Desireble_Segm_path)
    # after I export the segments for each segment I can perform k-means clustering to exrtract the seperate surfaces for each segment

#RE-SEGMENTATION
#Perform the segmentation process one more time, in the produced ...segments.ply file. Becasue now this file does not contain the noise that the initial file contained int he beginning
def ask_to_continue():
    while True:
        user_input = input("Do you want to continue with the Downsampling process? If Yes, press 1. If No, press 0: ")
        if user_input == '1':
            print("Continue Re-segmenting the PCM.")
            return True
        elif user_input == '0':
            print("Exiting the Re-segmenting process.")
            return False
        else:
            print("Invalid input. Please press 1 for Yes or 0 for No.")
#if ask_to_continue():
    #Multi-order RanSAC planar Segmentation
    #Take as an imput the first segmented PCM file
    ply_file_path_in2=  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segements.ply"
    pcd2 = read_pcd(ply_file_path_in2)
    print(pcd2)
    max_plane_index =44
    plane_to_pnt_dist = 0.02
    segment_models = {}
    segments = {}
    rest = pcd2
    Segments_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1//Segments/PLY_SDK_EU_Real.ply"
    for i in range(max_plane_index):
        colors = plt.get_cmap( "tab20")(i)
        segment_models[i], inliers = rest.segment_plane(distance_threshold = plane_to_pnt_dist, ransac_n = 3, num_iterations = 7000) 
        #we exclude inliers as indexes that can be interpreted only by open3d library's functions
        segments[i]= rest.select_by_index(inliers)
        segments[i].paint_uniform_color(list(colors[:3]))
        rest= rest.select_by_index(inliers, invert = True)
        directory, file_name_with_extension = os.path.split(Segments_path)
        file_name, extension = os.path.splitext(file_name_with_extension)
        additional_word = f"_Segment_{i+1}"
        # Modify the filename based on the condition
        new_file_name = file_name + additional_word
        # Recreate the output file paths
        Desirable_Segm_path_out = os.path.join(directory, new_file_name + extension)
        segment_points = np.asarray(segments[i].points)
        segment_colors = np.asarray(segments[i].colors)
        write_ply_from_pcd(segment_points, segment_colors*255, Desirable_Segm_path_out)

        print("pass", i,"/", max_plane_index, "done")
        #o3d.visualization.draw_geometries([inliers] )
    o3d.visualization.draw_geometries([segments[i]for i in range(max_plane_index)]+ [rest])
    #print(segments)
    #print(segment_models)
    txt_Segments_path = r"C:/Python Output/Real Data/TXT files//PLY_SDK_EU_Real_segements.txt"
    create_TXT(txt_Segments_path, segments)
    segm_points,segm_colors= read_txt(txt_Segments_path)
    Desireble_Segm_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segements_2.ply"
    write_ply_from_pcd(segm_points, segm_colors*255, Desireble_Segm_path)

#Interactive Segmentation Process
#In this process we apply the interactive Segmentation Algorithm.
#So we need a floor plan view of the newly segmented ply file
#We import the Py code that is repsonsible for denoising
#Denoising_py_file = r"C:/Python algorithm Thesis/Scripts Testing/Code 19_real dataset/3.9 SDK_real_Cube.py"
#subprocess.run(["python", Denoising_py_file])

#Translation of the Floor plan view ply file
Floor_plan_view_path=  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segements_2_Cube.ply"
Floor_plan_view_pcd = read_pcd(Floor_plan_view_path)
Origin_point = np.array([20,0, 0])
Floor_plan_Translated = Floor_plan_view_pcd.translate(Origin_point) 
points_tr = np.asarray(Floor_plan_Translated.points)
colors_tr = np.asarray(Floor_plan_Translated.colors)
Floor_plan_view_translated =  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segements_2_Cube_translated.ply"
write_ply_from_pcd(points_tr, colors_tr *255, Floor_plan_view_translated) #Scale colors back to the 0-255 range for PLY
# so now we have: 1) all the segments seperated in different ply files
                # 2) a ply file that represents the floor plan view, that is a horizontal section of the final segmented file
                # 3) a translated ply.file of the horizontal section
#We initialise again our ply files
#Floor_plan_view_pcd
#Floor_plan_Translated
#and the segments from the file
def Rewrite_segement_ply_files(Segment_Name, Seg_state, ply, output_directory):
    # Save the processed point cloud to a new file
    output_filename = f"{Segment_Name}_{Seg_state}.ply"
    output_filepath = os.path.join(output_directory, output_filename)
    o3d.io.write_point_cloud(output_filepath, ply)

all_point_clouds = [Floor_plan_view_pcd, Floor_plan_Translated]
input_Segement_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/'
output_Correct_Segment_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/Correct_Segments/'
output_Wrong_Segment_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/Wrong_Segments/'
output_Discarded_Segment_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/Discarded_Segments/'
Count_Correct_seg = 0
Count_Wrong_seg = 0
Count_Discarded_seg = 0
#Initialise an array in which we will store the the name of the segment if it is wrong or not an the real numbe of segments this segment should have been represented
#For example we may have a segment that should have been 2 segments, so we store in this array the number two. Then at the end we will sum all these numbers and we will run the Segmentation process again with the new number of segments.
Segment_array= np.array([["segment name", "State", "K-means clusters" ]], dtype=object)
merged_wrong_point_cloud = o3d.geometry.PointCloud()
for filename in os.listdir(input_Segement_directory):
    if filename.endswith('.ply'):
        # Load .ply file
        filepath = os.path.join(input_Segement_directory, filename)
        Segment_ply = o3d.io.read_point_cloud(filepath)
        Floor_plan_1 = o3d.io.read_point_cloud(Floor_plan_view_path)
        Floor_plan_2 = o3d.io.read_point_cloud(Floor_plan_view_translated)
        # Visualize the PLY files together
        o3d.visualization.draw_geometries([Floor_plan_1, Floor_plan_2, Segment_ply])
        print(Segment_ply)

        Boolean_i = float(input("Does this Segment represent correct the reality, If yes press 1, If no press 0, To discard this segments press 2"))
        Real_N_Seg = float(input("In how many Segments. this segment should be splited off?"))
        if Boolean_i == 1:
            Count_Correct_seg = Count_Correct_seg + 1
            Seg_state = "Correct"
            file_name, extension = os.path.splitext(filename)
            Rewrite_segement_ply_files(file_name, Seg_state, Segment_ply, output_Correct_Segment_directory)
        elif Boolean_i == 0:
            Count_Wrong_seg = Count_Wrong_seg + 1
            Seg_state = "Wrong"
            file_name, extension = os.path.splitext(filename)
            Rewrite_segement_ply_files(file_name, Seg_state, Segment_ply, output_Wrong_Segment_directory)
            merged_wrong_point_cloud += Segment_ply
        elif Boolean_i == 2:
            Count_Discarded_seg = Count_Discarded_seg +1
            Seg_state = "Discarded"
            file_name, extension = os.path.splitext(filename)
            Rewrite_segement_ply_files(file_name, Seg_state, Segment_ply, output_Discarded_Segment_directory)
        New_Seg_Object =  np.array([file_name, Seg_state, Real_N_Seg])
        Segment_array = np.vstack([Segment_array, New_Seg_Object])       
print(Count_Correct_seg)
print(Count_Wrong_seg)
print(Count_Discarded_seg)
merged_wrong_ply_points = np.asarray(merged_wrong_point_cloud.points)
merged_wrong_ply_colors = np.asarray(merged_wrong_point_cloud.colors)
Merged_wrong_seg =  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_Wrong_Segments.ply"
write_ply_from_pcd(merged_wrong_ply_points, merged_wrong_ply_colors*255, Merged_wrong_seg)

#Now we have created three folders that contain the Correct Sgments, the discarded Segments, and the Wrong Segments
#The wrong segments are the segments that contain one or more correct segments and some noise.
#First, step in dealing with those segments is to apply K-means clustering, to seperate the Segments in to different clusters. K-means Clustering is capable on finding groups of points that are close, and have similar geometrical characterisitcs.
#We use as an input the Wrong Segments that have been produced before.


        