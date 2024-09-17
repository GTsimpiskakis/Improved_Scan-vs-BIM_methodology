#B.2. Point cloud file Downsampling
#First this script executes ply file downsampling without any user's internvention

#B.3. Point cloud file Segmentation
#Second executes (denoised) file segmentation by applying RANSAC method
import numpy as np
import open3d as o3d
import matplotlib.pyplot as plt 
import os 

#pcd_file_path500k = r"C:/Users/dell/OneDrive/Υπολογιστής/CME/Season 2023-2024/Graduation Project/Thesis/Data/Test Data/Topoplan/Produced PCD from test code 1/Apartment.ply"
ply_file_path_in=  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_Cube.ply"
ply_file_path_out =  r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_DownSampled.ply"
pcd_file_path_out =  r"C:/Python Output/Real Data/PCD files/SDK_real_dataset_1/PLY_SDK_EU_Real_DownSampled.pcd"

def read_pcd(pcd_file_path):
    pcdfile = o3d.io.read_point_cloud(pcd_file_path)
    return pcdfile
pcd = read_pcd(ply_file_path_in)
print(pcd)

#Statistical outlier filter
nn = 8  #Number of nearest neighbours
std_multiplier = 3
filtered_pcd = pcd.remove_statistical_outlier(nn, std_multiplier)
outliers = pcd.select_by_index(filtered_pcd[1], invert = True)
outliers.paint_uniform_color([1,0,0])
filtered_pcd = filtered_pcd[0]
print(outliers)
o3d.visualization.draw_geometries([filtered_pcd, outliers])

#Voxelisation (Downsampling)
Voxel_size= 0.01
pcd_downsample = filtered_pcd.voxel_down_sample(voxel_size = Voxel_size)
o3d.visualization.draw_geometries([pcd_downsample])
# Save Open3D PointCloud object to PCD file
o3d.io.write_point_cloud(pcd_file_path_out, pcd_downsample, write_ascii=True)
print(pcd_downsample)

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

#Multi-order RanSAC planar Segmentation
max_plane_index =44
plane_to_pnt_dist = 0.015
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
Desireble_Segm_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/PLY_SDK_EU_Real_segements.ply"
write_ply_from_pcd(segm_points, segm_colors*255, Desireble_Segm_path)
# after I export the segments for each segment I can perform k-means clustering to exrtract the seperate surfaces for each segment