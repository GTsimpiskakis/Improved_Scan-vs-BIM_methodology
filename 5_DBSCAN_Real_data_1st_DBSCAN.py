import numpy as np
import open3d as o3d
import matplotlib.pyplot as plt 
import os

def read_pcd(pcd_file_path):
    pcdfile = o3d.io.read_point_cloud(pcd_file_path)
    return pcdfile

def read_pcd2(pcd_file_path):
    pcd = o3d.io.read_point_cloud(pcd_file_path)
    points = np.asarray(pcd.points)
    colors = np.asarray(pcd.colors)
    return points, colors
                                                      
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
pcd_file_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/Wrong_Segments/PLY_SDK_EU_Real_Segment_12_Wrong.ply"
Input_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/DBSCAN_clusters'

pcd = read_pcd(pcd_file_path)
points, colors= read_pcd2(pcd_file_path)#read the downsampled pcd
#write_ply_from_pcd(points, colors *255, ply_file_path_out) #Scale colors back to the 0-255 range for PLY
print("Conversion completed.")# conversion from downsampled pcd to ply completed


directory, file_name_with_extension = os.path.split(pcd_file_path)
filename, extension = os.path.splitext(file_name_with_extension)
# Split the file name by underscore '_'
file_name_parts = filename.split('_')
# The value "11" is the 5th element in the list (index 4)
segment_number = file_name_parts[5]
folder_filename = f"{file_name_parts[0]}_{file_name_parts[1]}_{file_name_parts[2]}_{file_name_parts[3]}_{file_name_parts[4]}_{file_name_parts[5]}"
# Define the new output directory and create a new folder in this for storing the PLY files, the new folder's name is filename
Output_directory = os.path.join(Input_directory, folder_filename)
# Use os.makedirs to create the new folder
os.makedirs(Output_directory, exist_ok=True)
#Create two more folders which will be used to classify manualy the produced segments from DBSCAN
Output_directory_2 = os.path.join(Output_directory, "Correct_Segments")
os.makedirs(Output_directory_2, exist_ok=True)
Output_directory_3 = os.path.join(Output_directory, "Discarded_Segments")
os.makedirs(Output_directory_3, exist_ok=True)
#Apply DBSCAN algoirthm
with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as cm:
    labels = np.array(pcd.cluster_dbscan(eps=0.02, min_points=5, print_progress=True))

max_label = labels.max()
print(f"point cloud has {max_label + 1} clusters")
colors = plt.get_cmap("tab20")(labels / (max_label if max_label > 0 else 1))
colors[labels < 0] = 0
pcd.colors = o3d.utility.Vector3dVector(colors[:, :3])
# Save each cluster to a separate PLY file
for cluster_label in range(max_label + 1):
    cluster_indices = np.where(labels == cluster_label)[0]
    cluster_pcd = pcd.select_by_index(cluster_indices)
    output_filename = f"{file_name_parts[0]}_{file_name_parts[1]}_{file_name_parts[2]}_{file_name_parts[3]}_{file_name_parts[4]}_{file_name_parts[5]}_{cluster_label}.ply"
    cluster_pcd_path = os.path.join(Output_directory, output_filename)
    o3d.io.write_point_cloud(cluster_pcd_path, cluster_pcd)

# We visualize the entire point cloud with clusters to manually inspect if the RANSAC segment has been split in a way that the desirable segments have been produced
o3d.visualization.draw_geometries([pcd])
                                  
